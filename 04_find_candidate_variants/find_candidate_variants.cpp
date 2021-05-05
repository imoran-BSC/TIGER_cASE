/*
 *  FIND CANDIDATE VARIANTS
 *
 *  Reads an Allelic_Imba_trimmed.txt file and the full genotype
 *  of the given chr, and for each reporter SNP, calculates
 *  the Zscores from each genotyped SNP (separating Het and Hom
 *  samples) within a TAD+leeway.
 *
 *  Writen by Ignasi Moran
 */

#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <map>
#include <set>

#include "../utils/class_fstream.h"
#include "../utils/allelic_imba_common.h"  // Contains the functions to calculate imba, covg criteria
// ../utils/allelic_imba_common.cpp  also needs to be included in the project or compile command
// since it contains the implementations of the functions defined in the header.


using namespace std;

bool debug = 0,  // DEBUG FLAG
  no_limits = 0; // Remove ALL pair limits, to fetch otherwise discarded "interesting" pairs
const string testchr = "chr21";

const string target = "";
// To debug specific rsids if necessary

const long int tad_leeway = 200000;  // Expand the search this much around TAD borders
// MIGHT BE WORTH TO EXPAND 1 FULL TAD INSTEAD
const double min_r2_limit = 0.8,  // Min r2 value for value to be considered
  min_maf_value = 0.005,
  hom_pval_limit = 0.05,  // Max allowed p-value for homozygous samples
  hwe_pval_limit = 0.001;  // Limit to reject a snp for HWE disequilibrium,
const string cromosomes[24] =
  {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
   "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
   "chr20","chr21","chr22","chrX","chrY"};

/*
 *   CLASSES
 */

double funct_calc_pval( double zval, vector< pair<double,double> >& distr, int tails );  // Declaration
double funct_calc_zvalue_np( vector<double*>& pv, vector<double>& ratios,
  vector<int*>& read_counts, int with_direction=1 );

class class_region
{  // Contains a genomimc bed region and its type
public:
  long int ini, fin;
  string chr, type;
  list<string> genesymb_list;

  class_region() : ini(0), fin(0), chr(""), type("") {}

  void read( class_fstream& in )
    { in.file >> chr >> ini >> fin; }

  void read_tads( class_fstream& in )
    {  // Reads TAD file where type=0 is intraTAD space, and type=1 is a TAD
    in.file >> chr >> ini >> fin >> type;
    // Now remove tad leeway to nonTADs not have pairs both in and out of TADs
    if( type == "0" ) {  // Non TADs
      ini += 2*tad_leeway+1; fin -= 2*tad_leeway+1;  // This will get +tad_leeway later in analysis
      if( fin+tad_leeway < ini-tad_leeway ) { ini = -1; fin = -1; } }
    }

  bool operator< ( const class_region& r ) const
    { if(ini != r.ini) return ini < r.ini; else return fin < r.fin; }
};

class class_genotype_phased
{  // Copied from class genotype, modified to deal with phasing
public:
  long int pose;
  string chr, rsname, ref, alt, reg_overlap,
    gwas_ld_rsnames, gwas_ld_r2vals, gwas_ld_type, eqtl_study;
  vector<string> genoty_vec;
  list<string> genesymb_list;

  class_genotype_phased() : pose(0), chr(""), rsname(""), ref(""), alt(""), reg_overlap(""),
    gwas_ld_rsnames(""), gwas_ld_r2vals(""), gwas_ld_type(""), eqtl_study("") {};

  void read( class_fstream& in )
    {  // Llegeix la info de la entrada i simplifica el genotipatge
    string dummy="", gty="", prob="";
    vector<string> nt_vec( num_genoty_samples );
    genoty_vec.clear(); genoty_vec.resize( num_genoty_samples );

    in.file >> chr >> pose >> rsname >> ref >> alt >> dummy >> dummy >> dummy >> dummy;
    chr = "chr"+chr;

    for( unsigned k=0; k<num_genoty_samples; k++ )
      {  // Converteix el genotipatge de NTs a Ref/Het/Alt/NAN
      in.file >> nt_vec[k];
      gty=nt_vec[k].substr( 0, nt_vec[k].find(":") );
      prob=nt_vec[k].substr( nt_vec[k].find(":")+1 );

      genoty_vec[k] = resolve_genotypes( gty, prob );
      }
    }

  string resolve_genotypes( string& gty, string& p )
    {  // Given a genotype and its probability, resolves the phased genotype
    string res="";
    // Do NOT take "genomic dosage" as a quality measure (prob), since it performs way worse if we do
    if( gty == "0|0" || gty == "0/0" ) res = "Ref";
    else if ( gty == "0/1" ) res = "Het";
    else if ( gty == "0|1" ) res = "Het1";  // Distinguishes the two phased heterozygous
    else if ( gty == "1|0" ) res = "Het2";
    else if( gty == "1|1" || gty == "1/1" ) res = "Alt";
    else res = "NAN";

    return res;
    }

  bool operator< ( const class_genotype& g ) const
    { return pose < g.pose; }
};

class class_snp_pair
{  // Contains all info and functions of a reporter-cisreg pair
public:
  unsigned num_good_hets, num_good_homs, num_het_signif, num_hom_signif;
  long int dist;
  double r2val, zvalue, zpvalue, repz, jointar, homzvalue, homzpvalue,
    abf, pprob;
  string credset;

  vector<int> phase_vec;   // 0 for same phase, 1 for opposite
  vector<int*> hets_readc_pvec, homs_readc_pvec, hets_refaltc_pvec;
  vector<double*> hets_imba_pvec, hets_pval_pvec, homs_imba_pvec, homs_pval_pvec;  // Het and Hom GENOTY
  vector<string*> hets_sample_pvec, homs_sample_pvec;

  class_snp* snpp;  // Reporter SNP
  class_genotype_phased* genotyp;  // Candidate causal SNP
  class_genotype_phased genoty;    // Actual genotype data copy
  class_region* tadp;

  class_snp_pair() : num_good_hets(0), num_good_homs(0), num_het_signif(0),
    num_hom_signif(0), dist(0), r2val(0), zvalue(0), zpvalue(0), repz(0),
    jointar(0), homzvalue(0), homzpvalue(0), abf(0), pprob(0),
    credset(""), genotyp(0) {}

  bool operator< (const class_snp_pair& s ) const
    {  // Sorts the SNP pairs by pvalue then zscore then r2
    if( zpvalue != s.zpvalue ) return zpvalue < s.zpvalue;   // Sorts by increasing zpvalue
    else if( fabs(zvalue) != fabs(s.zvalue) ) return fabs(zvalue) > fabs(s.zvalue);  // Breaks the ties by higher Zscore
    else return r2val > s.r2val;  // And then by r2 value
    }

  void clean()
    {  // Clears the vectors
    genotyp=0;
    hets_imba_pvec.clear(); hets_pval_pvec.clear(); homs_imba_pvec.clear(); homs_pval_pvec.clear();
    hets_sample_pvec.clear(); homs_sample_pvec.clear();
    phase_vec.clear(); hets_readc_pvec.clear(); homs_readc_pvec.clear();
    hets_refaltc_pvec.clear();
    }

  void calc_zvalues( vector< pair<double,double> >& rndzval_distr )
    {  // Calculates the zvalues and zpvalues for het and hom cisreg populations
    vector<double> hets_imba_phcand_vec( phase_vec.size() );

    for( unsigned u=0; u<phase_vec.size(); u++ )
      {  // To adress the phasing
      if( phase_vec[u] == 0 ) hets_imba_phcand_vec[u] = *(hets_imba_pvec[u]);
      else if( phase_vec[u] == 1 ) hets_imba_phcand_vec[u] = 1 - *(hets_imba_pvec[u]);
      else { cout << "What? " << phase_vec[u] << " found!" << endl; cin.get(); }
      }

    zvalue = funct_calc_zvalue_np( hets_pval_pvec, hets_imba_phcand_vec, hets_readc_pvec );
    homzvalue = funct_calc_zvalue( homs_pval_pvec, homs_imba_pvec, homs_readc_pvec );

    zpvalue = funct_calc_pval( zvalue, rndzval_distr, 2 );
    homzpvalue = funct_calc_pval( homzvalue, rndzval_distr, 2 );
    }

  void make_pair( class_genotype_phased& g, class_snp& s,
      class_region* tp )
    {  // Copy the genotype info since it is NOT stored in RAM
    map< string, double >::iterator benj_it;
    num_het_signif=0; num_hom_signif=0;

    // Calculate the missing output values here
    genotyp = &g; snpp = &s; tadp = tp;
    dist = genotyp->pose - snpp->pose;  // Distance from the reporter SNP

    // Count the number of consistent and inconsistent samples
    num_good_hets = hets_pval_pvec.size(); num_good_homs = homs_pval_pvec.size();
    for(unsigned u=0; u<num_good_hets; u++)
      if( *(hets_pval_pvec[u]) < binom_pval_threshold )
        num_het_signif++;
    for(unsigned u=0; u<num_good_homs; u++)
      if( *(homs_pval_pvec[u]) < binom_pval_threshold )
        num_hom_signif++;  // Inconsist as in significant HOM candidates
    }

  void copy_genoty()  // Expensive, only done for pairs pushed to the vec
    { genoty = *genotyp; genotyp = 0; }  // Need to copy since genoty will expire due to looping

  bool self_pair()
    {  // pair with itself (Reporter)? - Only works after copy_genoty!!
    return genoty.rsname == snpp->rsname;
    }

  bool candidate_Z()
    {  // Check if candZ >= repZ, only works after having assigned repz!!
    return fabs( zvalue ) >= fabs( repz );
    }

  string string_nsignif_samples()
    {  // Returns number of ref/alt biased samples in a string
    unsigned num_ref_bias=0, num_alt_bias=0;
    stringstream ss;

    for( unsigned k=0; k<num_good_hets; k++ )
      if( *(hets_pval_pvec[k]) <= binom_pval_threshold )
        { if( *(hets_imba_pvec[k]) > 0.5 ) num_ref_bias++;
          else num_alt_bias++; }

    ss << num_ref_bias << "|" << num_alt_bias;
    return ss.str();
    }

  void write_reporter( class_fstream& out )
    {  // Write the reporter-centric Zscore info
    out.file << genoty.rsname << "\t" << zvalue << "\t" << zpvalue << "\t"
      << jointar << "\t" << snpp->median_coverage << "\t"
      << snpp->chr << "\t" << snpp->pose << "\t";
    out.file << string_nsignif_samples() << "\t";

    if( num_good_hets > 0 )
      {  // Het info
      for( unsigned k=0; k<num_good_hets-1; k++ ) out.file << *(hets_imba_pvec[k]) << ",";
      out.file << *( hets_imba_pvec.back() ) << "\t";
      for( unsigned k=0; k<num_good_hets-1; k++ ) out.file << *(hets_pval_pvec[k]) << ",";
      out.file << *( hets_pval_pvec.back() ) << "\t";
      for( unsigned k=0; k<num_good_hets-1; k++ ) out.file << *(hets_sample_pvec[k]) << ",";
      out.file << *( hets_sample_pvec.back() ) << "\t";
      // Added Refc/Alt raw read counts as well
      out.file << *(hets_refaltc_pvec[0]) << "|" << *(hets_refaltc_pvec[1]);
      for( unsigned k=1; k<num_good_hets; k++ )
        out.file << "," << *(hets_refaltc_pvec[2*k]) << "|"
          << *(hets_refaltc_pvec[2*k+1]);
      }
    else out.file << "NA\tNA\tNA\tnone\tNA";
    out.file << "\n";
    }

  void writeout( class_fstream& out )
    {  // Writes out all the relevant info of the pair
    out.file << genoty.rsname << "\t" << snpp->rsname << "\t" << r2val << "\t"
      << snpp->median_coverage << "\t"
    // Then all relevant info about the reporter SNP/Z-score
      << zvalue << "\t" << zpvalue << "\t" << repz << "\t"
      << abf << "\t" << pprob << "\t" << credset << "\t"
      << jointar << "\t" << snpp->chr << "\t" << snpp->pose << "\t"
      << genoty.pose << "\t" << dist << "\t";
    // TAD info
    if( tadp->type == "1" )
      out.file << tadp->ini << "\t" << tadp->fin << "\t" << tadp->type << "\t";
    else out.file << tadp->ini-tad_leeway << "\t" << tadp->fin+tad_leeway << "\t"
        << tadp->type << "\t";
    // Detailed data and sample-specific values
    out.file << num_het_signif << "\t" << num_hom_signif << "\t"
      << homzvalue << "\t" << homzpvalue << "\t";
    if( num_good_hets > 0 )
      {  // Het info
      // Prints the UNPHASED AR values
      for( unsigned k=0; k<num_good_hets-1; k++ ) out.file << *(hets_imba_pvec[k]) << ",";
      out.file << *( hets_imba_pvec.back() ) << "\t";
      for( unsigned k=0; k<num_good_hets-1; k++ ) out.file << *(hets_pval_pvec[k]) << ",";
      out.file << *( hets_pval_pvec.back() ) << "\t";
      for( unsigned k=0; k<num_good_hets-1; k++ ) out.file << *(hets_sample_pvec[k]) << ",";
      out.file << *( hets_sample_pvec.back() ) << "\t";
      // Added Refc/Alt raw read counts as well
      out.file << *(hets_refaltc_pvec[0]) << "|" << *(hets_refaltc_pvec[1]);
      for( unsigned k=1; k<num_good_hets; k++ )
        out.file << "," << *(hets_refaltc_pvec[2*k]) << "|"
          << *(hets_refaltc_pvec[2*k+1]);
      out.file << "\t";
      }
    else out.file << "NA\tNA\tNA\tnone\tNA\t";
    if( num_good_homs > 0 )
      {  // Hom info
      for( unsigned k=0; k<num_good_homs-1; k++ )
        out.file << *(homs_imba_pvec[k]) << ",";
      out.file << *( homs_imba_pvec.back() ) << "\t";
      for( unsigned k=0; k<num_good_homs-1; k++ ) out.file << *(homs_pval_pvec[k]) << ",";
      out.file << *( homs_pval_pvec.back() ) << "\t";
      for( unsigned k=0; k<num_good_homs-1; k++ ) out.file << *(homs_sample_pvec[k]) << ",";
      out.file << *( homs_sample_pvec.back() ) << "\t";
      }
    else out.file << "NA\tNA\tnone";
    out.file << "\n";
    }
};

/*
 *   TERTIARY FUNCTIONS
 */

bool funct_sort_rsids( const class_snp_pair& p1, const class_snp_pair& p2 )
{ if( p1.snpp->rsname != p2.snpp->rsname ) return p1.snpp->rsname.compare( p2.snpp->rsname ) > 0;
  else return p1.genoty.rsname.compare( p2.genoty.rsname ) > 0; }  // Sorts high to low coverage

bool funct_unique_rsids( const class_snp_pair& p1, const class_snp_pair& p2 )
{ return p1.snpp->rsname == p2.snpp->rsname && p1.genoty.rsname == p2.genoty.rsname; }

bool funct_sort_abfp( const class_snp_pair* p1, const class_snp_pair* p2 )
{ return p1->abf > p2->abf; }  // High to low ABF pair sorting

void funct_load_tads( class_fstream& in, vector< vector< class_region > >& templ_vec )
{  // Reads the tad region using type ISTAD, since we now have NEGATIVE TADs w type 0, TADs are 1
unsigned chrn=0;
string chrt="";

class_region templ;

templ.read_tads( in );
while( !in.file.eof() )
  {  // Reads all (UCSC, RefSeq, Ensembl) and lncRNA exons
  if( chrt != templ.chr )
    { chrt=templ.chr; for(chrn=0; chrn<22; chrn++) if( cromosomes[chrn] == chrt ) break; }

  if( chrn<22 && templ.fin > 0 ) templ_vec[chrn].push_back( templ );

  templ.read_tads( in );
  }
}

void funct_load_r2( class_fstream& in, map<string, double>& r2_map )
{  // Loads the r2 map, encoded as rsid1+rsid2 -> r2 value !!! (was - before!!)
double r2=0;
string dummy, rs1, rs2;

in.file >> dummy >> dummy >> rs1 >> dummy >> dummy >> rs2 >> r2;
while( !in.file.eof() )
  {
  if( r2 >= min_r2_limit )  // Only load to memory pairs over r2 limit
    r2_map[ rs1 + "+" + rs2 ] = r2;

  in.file >> dummy >> dummy >> rs1 >> dummy >> dummy >> rs2 >> r2;
  }
}

void funct_load_maf( class_fstream& in, map<long int, double>& maf_map )
{  // Loads the MAF map
size_t p1=0, p2=0;
long int pose=0;
double maf=0;
string dummy, rsid;

in.file >> dummy >> rsid >> dummy >> dummy >> maf >> dummy;
while( !in.file.eof() )
  {
  p1 = rsid.find(":");
  p2 = rsid.find("_");
  pose = atof( rsid.substr( p1+1, p2-p1+1 ).c_str() );

  if( maf >= min_maf_value )
    maf_map[ pose ] = maf;

  in.file >> dummy >> rsid >> dummy >> dummy >> maf >> dummy;
  }
}

void funct_load_hwe( class_fstream& in, set<string>& nonhw_rsid_set )
{  // Loads the nonHWE rsids set
double pval=0;
string dummy, rsid;

in.file >> dummy >> rsid >> dummy >> dummy >> dummy >> dummy
  >> dummy >> dummy >> pval;
while( !in.file.eof() )
  {
  if( pval <= hwe_pval_limit )
    nonhw_rsid_set.insert( rsid );

  in.file >> dummy >> rsid >> dummy >> dummy >> dummy >> dummy
  >> dummy >> dummy >> pval;
  }
}

double funct_calc_zvalue_np( vector<double*>& pv, vector<double>& ratios,
  vector<int*>& read_counts, int with_direction )  // with_direction defaults to 1 in declaration
{  // Carbon copy but with vector<double>& ratios instead of a double*
unsigned good=0, vecsize = pv.size();
double d=0, up=0, down=0;

vector<double> z(vecsize);

for( unsigned i=0; i<vecsize; i++ )
  {
  if( *(pv[i]) <= 0 || ratios[i] <=0 || *(read_counts[i]) <= 0 ) continue;  // Skip this one - should not happen

  z[i] = funct_qnorm( *(pv[i])/2.0 );  // Fetches the right Z-scores from the R table, z=qnorm(1-p/2)

  d=1;
  if( with_direction && ratios[i] < 0.5 ) d = -1;
  if( ratios[i] == 0.5 ) d = 0;  // Remove the exactly 0.5 ARs
  up += z[i]*sqrt( *(read_counts[i]) )*d;  // as per paper Z = sum(z_i * w_i * d_i )/sqrt( sum(w_i^2) )
  down += *(read_counts[i]);
  good++;
  }

if( down > 0 && good >= min_num_hets ) return up/sqrt( down );
else return -100;
}


/*
 *   SECONDARY FUNCTIONS
 */

void blob_of_text()
{  // A blob of text for when the wrong arguments are passed
cout << "\nThis program reads an Allelic_Imba_trimmed.txt file "
  << "calculates the z-scores of all its SNPs, and reports their significance.\n\n"
  << "It requires 3 arguments:\n"
  << " 1) The general path to where the Allelic_Imba_trimmed.txt "
  << "(output from allelic_imba_filter_and_unify) and "
  << "Allelic_Imba_rnd_Zscore.txt files are located, ending in /. "
  << "Ie /project/jferrer/imperial/SNP_Analysis/\n"
  << " 2) The general path to where the gene (Genes_UCSC_RefSeq_Ensembl_lncRNAs_hg19.bed) "
  << "and exon (Exons_UCSC_RefSeq_Ensembl_lncRNAs_hg19.bed) annotation files are located"
  << "plus the geneSymbol conversion files (Convert_UCSC2GSymb_hg19.txt, Convert_RefSeq2GSymb_hg19.txt)"
  << "and Convert_Ensembl2GSymb_hg19.txt), plus the Regulome and GWAS folders, "
  << "which should contain the regulome files (C{1..5}_sites.bed and clustered_C3_sites.bed) "
  << "and the gwas files (GWAS_FG.bed, GWAS_T2G.bed). Ie /project/jferrer/imperial/Genomic_Info/\n"
  << " 3) The path where the output files will be located, ending in /. "
  << "Ie /project/jferrer/imperial/Fastq_Data/Results/\n"
  << "\nFor example, a valid call would be\n"
  << "$ ./allelic_imba_zscore /project/jferrer/imperial/SNP_Analysis/ "
  << "/project/jferrer/imperial/Genomic_Info/ "
  << "/project/jferrer/imperial/Fastq_Data/Results/\n"
  << "\nIf any of the above expectations are not met, the program will most probably "
  << "simply crash without much explanation or get stuck in an infinite loop.\n"
  << "\nThere is one output file:\n"
  << " - Allelic_Imba_Zscores.txt , which contains, for all SNPs in the input file, "
  << "the z-score, its p-value and Benjamini-Hochberg significance, plus information on which samples "
  << "was the SNP imbalanced, which genes intersected, which GWAS regions, etc.\n\n"
  << "This software is provided without any warranty whatsoever, so use it at your own risk.\n";
}

double funct_calc_pval( double val, vector< pair<double,double> >& distr, int tails )
{  // Given the val and the null distribution (and 1 or 2 tails), return the pvalue
unsigned distr_size=distr.size();
double pval=0, cutoff=0, interval = distr[2].first-distr[1].first;

if( val == 0 || val <= -100 ) return 1;

if( tails == 2 )
  {  // Two tailed, so we search for abs(val) and -1*abs(val) and add up both values
  if(val >= 0) {  // Search for the cutoff from the midpoint of the distr
    for( unsigned w=distr_size/(double)2-1; w<distr_size; w++ )
      if( val < distr[w].first )  // Search only positive space
        { cutoff = distr[w-1].second; break; } }  // Get the cutoff value
  else for( unsigned w=distr_size/(double)2+1; w>0; w-- )
    if( val > distr[w].first )  // Search only negative space
      { cutoff = distr[w+1].second; break; }  // Get the cutoff value
  for( unsigned w=0; w<distr_size; w++ )
    if( distr[w].second <= cutoff )
      pval += distr[w].second;  // Add up all values less than that (2-tailed)
  }
else if( tails == 1)
  {  // One tailed for Abs distribution
  for( unsigned w=1; w<distr_size; w++ )
    if( val <= distr[w].first )
      pval += distr[w-1].second;  // Adds all 1tail values after it (except last zero..)
  }
else { cout << "Can't compute a pvalue with " << tails << " tails!" << endl; cin.get(); }

pval = pval*interval;

if( pval < min_pval ) pval = min_pval;  // Lowest possible value
return pval;
}

void funct_read_rndzscore_and_sort( class_fstream& rnd_zval_in,
  vector< pair<double,double> >& rndzval_distr )
{  // Loads the control distribution, sorting the values to guarantee smooth continuity
string dummy="";
double dumz=0, dump=0;
vector<double> zval_neg, zval_pos, pval_neg, pval_pos;

rnd_zval_in.file >> dumz >> dummy >> dump;
while( !rnd_zval_in.file.eof() )
  {
  if( dumz < 0 ) { zval_neg.push_back(dumz); pval_neg.push_back(dump); }
  else { zval_pos.push_back(dumz); pval_pos.push_back(dump); }
  rnd_zval_in.file >> dumz >> dummy >> dump;
  }
sort( pval_neg.begin(), pval_neg.end() );
sort( pval_pos.rbegin(), pval_pos.rend() );  // Reverse iterators to sort in reverse
for( unsigned k=0; k<pval_neg.size(); k++ )
  rndzval_distr.push_back( pair< double, double >( zval_neg[k], pval_neg[k] ) );
for( unsigned k=0; k<pval_pos.size(); k++ )
  rndzval_distr.push_back( pair< double, double >( zval_pos[k], pval_pos[k] ) );
}

void funct_load_tad_snps( class_region& tad, list<class_snp>& repsnp_list,
  list<class_snp*>& tad_snpps )  // 'Fast' copy of iterator
{  // Loads the SNPs overlapping the TAD in the tad_snps list - added LEEWAY to the TAD boundaries
list<class_snp>::iterator snp_it;
tad_snpps.clear();

// Could be made faster by keeping the pointer to the first snp in the previous TAD but who cares
for( snp_it = repsnp_list.begin(); snp_it != repsnp_list.end() && (*snp_it).pose < tad.fin + tad_leeway; ++snp_it )
  if( (*snp_it).pose >= tad.ini - tad_leeway ) tad_snpps.push_back( &(*snp_it) );
}

void funct_search_pairs_and_zscores( vector< class_genotype_phased >& phcand_vec,
  list<class_snp*>& repsnp_list, vector<class_region>::iterator& tad_it,
  vector< int >& genoty2rnaseq_indx_vec, map<string, double>& r2_map,
  vector< pair<double,double> >& rndzval_distr,
  vector<string>& rnaseq_vec,
  vector<class_snp_pair>& rep_vec,
  map< string, double >& rep_zscore_map,
  vector<class_snp_pair>& pair_vec )
{  // Computes all Zscores of reporter SNPs within the TAD for the given genotpyed gene
unsigned jointp1c=0, jointp2c=0, phase_size = phcand_vec.size();
string dummy="";
map<string, double>::iterator r2_it, map_it;
set<string>::iterator set_it1, set_it2;
list<class_snp*>::iterator snp_it;
pair<double,double> pr;
class_genotype_phased *rep_phasep=0, *cas_phasep=0;

class_snp_pair spair;

for( snp_it = repsnp_list.begin(); snp_it != repsnp_list.end(); ++snp_it )
  {  // Iterates all reporter SNPs in the TAD
  dummy = "not_found";
  // First find the genotype for the reporter SNP

  for( unsigned uu=0; uu<phase_size; uu++ )
    if( phcand_vec[uu].pose >= (*tad_it).ini - tad_leeway &&
        phcand_vec[uu].pose <  (*tad_it).fin + tad_leeway &&
        (*snp_it)->rsname == phcand_vec[uu].rsname )
      { rep_phasep = &(phcand_vec[uu]); dummy = "found"; break; }
  if( dummy == "not_found" ) continue;

  for( unsigned uu=0; uu<phase_size; uu++ )
    if( phcand_vec[uu].pose >= (*tad_it).ini - tad_leeway &&
        phcand_vec[uu].pose <  (*tad_it).fin + tad_leeway )
    {  // Iterates all potential causal SNPs and makes pairs if they fit
    cas_phasep = &(phcand_vec[uu]);

    // Find which samples are het and their phasing for the causal SNP
    // Iterate all reporter SNPs, and calculate their het-only samples' Zscores
    spair.clean(); jointp1c=0, jointp2c=0;  // parent1 and parent2 read counts (phased)

    for( unsigned k=0; k<num_rnaseq_samples; k++ )
      if( genoty2rnaseq_indx_vec[k] != -1 && funct_coverage_test( *(*snp_it), k ) )
        {  // If it abides current read depth standard, get the corresponding Het values
        if( ( cas_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] == "Het"  ||
              cas_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] == "Het1" ||
              cas_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] == "Het2" ) &&
            ( rep_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] == "Het"  ||
              rep_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] == "Het1" ||
              rep_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] == "Het2" ) )
          {
          if( cas_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] ==
                rep_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] )
            { spair.phase_vec.push_back( 0 );  // Same phase
            if( (*snp_it)->sample_vec[ k ].pval < binom_pval_threshold ) {
              jointp1c += (*snp_it)->sample_vec[ k ].refc;
              jointp2c += (*snp_it)->sample_vec[ k ].altc; } }
          else  // If no phased, ASSUME all are in-phase !
            { spair.phase_vec.push_back( 1 );  // Opposite phase
            if( (*snp_it)->sample_vec[ k ].pval < binom_pval_threshold ) {
              jointp1c += (*snp_it)->sample_vec[ k ].altc;
              jointp2c += (*snp_it)->sample_vec[ k ].refc; } }

          spair.hets_imba_pvec.push_back( &(*snp_it)->sample_vec[ k ].imba );
          spair.hets_pval_pvec.push_back( &(*snp_it)->sample_vec[ k ].pval );
          spair.hets_readc_pvec.push_back( &(*snp_it)->sample_vec[ k ].totalcount );
          spair.hets_sample_pvec.push_back( &rnaseq_vec[ k ] );

          // We also make a refc+altc vector (2x as long) - redundant w readc but whatever
          spair.hets_refaltc_pvec.push_back( &(*snp_it)->sample_vec[ k ].refc );
          spair.hets_refaltc_pvec.push_back( &(*snp_it)->sample_vec[ k ].altc );
          }
        else if( cas_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] == "Ref" ||
            cas_phasep->genoty_vec[ genoty2rnaseq_indx_vec[k] ] == "Alt" )
          {  // If it abides current read depth standard, get the corresponding Hom values
          spair.homs_imba_pvec.push_back( &(*snp_it)->sample_vec[ k ].imba );
          spair.homs_pval_pvec.push_back( &(*snp_it)->sample_vec[ k ].pval );
          spair.homs_readc_pvec.push_back( &(*snp_it)->sample_vec[ k ].totalcount );
          spair.homs_sample_pvec.push_back( &rnaseq_vec[ k ] );
          }
        }

    spair.make_pair( *cas_phasep, *(*snp_it), &(*tad_it) );  // Copy all the missing info

    if( no_limits ||
        (*snp_it)->rsname == cas_phasep->rsname ||  // Its a reporter self-pair
        ( spair.num_hom_signif <= min_num_hets &&  // A good candidate pair
          spair.num_het_signif >= min_num_hets ) )
      {  // If not enough matching hets, we go to next causal
      // Search for the r2 value of the pair; if not in the map (lower than 0.8/rsid not present) then -1.
      spair.r2val= 0;
      r2_it = r2_map.find( cas_phasep->rsname + "+" + (*snp_it)->rsname );
      if( r2_it != r2_map.end() ) spair.r2val = (*r2_it).second;
      else { r2_it = r2_map.find( (*snp_it)->rsname + "+" + cas_phasep->rsname );  // Search the reverse
        if( r2_it != r2_map.end() ) spair.r2val = (*r2_it).second; }
      // We get the r2 value if returned by PLINK, otherwise zero
      // There was a -1 option before for non-calculated SNPs, now removed

      if( jointp1c+jointp2c>0 ) spair.jointar = jointp1c/(double)(jointp1c+jointp2c);
      else spair.jointar = -1;

      spair.copy_genoty();  // Commit the genoty data to RAM
      spair.calc_zvalues( rndzval_distr );

      // Add the reporter to its vector
      if( spair.zvalue > -100 )
        {
        if( spair.self_pair() )  // self pair now works w copied genoty
          {
          map_it = rep_zscore_map.find( spair.snpp->rsname );
          if( map_it == rep_zscore_map.end() )
            {  // First time we see the reporter
            rep_zscore_map[ spair.snpp->rsname ] = spair.zvalue;
            rep_vec.push_back( spair );  // Push the reporter self-pair
            }
          }

        // CONDITIONS FOR CONSIDERING THE PAIR ARE **HERE**
        // now ALSO includes self-pairs!!
        if( no_limits ||  // And only the worthy candidate pairs
            ( spair.homzvalue == -100 ||
              spair.homzpvalue > hom_pval_limit ) )
          pair_vec.push_back( spair );
        }
      }
    }
  }
}

void funct_write_reporters( vector<class_snp_pair>& rep_vec,
  class_fstream& out )
{  // Makes a reporter zscore map, and writes out all consistent reporters
unsigned goodc=0, badc=0;
vector<class_snp_pair>::iterator vec_it;

sort( rep_vec.begin(), rep_vec.end() );  // Sort by ascending Z-score pvals

for( vec_it=rep_vec.begin(); vec_it!=rep_vec.end(); ++vec_it )
  {
  if( (*vec_it).num_good_hets >= min_num_hets )
    {  // Writeout
    (*vec_it).write_reporter(out);
    goodc++;
    }
  else badc++;
  }

cout << "Wrote " << goodc << " good and not " << badc
  << " reporters out of " << rep_vec.size() << "\n";
}

void funct_calc_abf( vector<class_snp_pair*>& pair_pvec,
  map<long int, double>& maf_map, map< string, double >& rep_zscore_map )
{  // Given the set of SNP pairs, calculate the ABF -- Needs the SNP MAF!!
unsigned pvec_size = pair_pvec.size();
int casen=0, controln=0;  // Used # samples to compute the Zscore
double maf=0, revmaf=0, v=0, r=0, w=0.04,  // As per Silvia / SNPTEST
  credset_limit = 0.95,  // By default calculate 95% credsets
  cumul_abf=0, cumul_pprob=0, temp_Z=0,
  temp_abf=0, max_abf = 1E299;  // Attempt to remove Inf ABF vals?
map<long int, double>::iterator map_it;
map< string, double >::iterator rep_it;

for(unsigned u=0; u < pvec_size; u++)
  {
  map_it = maf_map.find( (pair_pvec[u])->genoty.pose );
  if( map_it != maf_map.end() ) maf = (*map_it).second;
  else maf = -1;

  // Fetch reporter Z for the pair, and store it in spair
  rep_it = rep_zscore_map.find( (pair_pvec[u])->snpp->rsname );
  if( rep_it != rep_zscore_map.end() )
    (pair_pvec[u])->repz = rep_zscore_map[ (pair_pvec[u])->snpp->rsname ];
  else cout << "Got NO reporter Z for " << (pair_pvec[u])->snpp->rsname << "\n";

  // We are using Z DELTA wrt reporter SNP, which just removes a constant
  // for all snp pairs sharing a reporter - keeps ABFs more reasonable
  if( (pair_pvec[u])->candidate_Z() )
    {
    temp_Z = fabs( (pair_pvec[u])->zvalue ) - fabs( (pair_pvec[u])->repz );

    if( maf != -1 )
      {
      revmaf = 1-maf;
      casen = num_rnaseq_samples;
      controln = casen;  // Just all equal 105
      // Since num_good_homs is zero in many cases which only have het-het and hom-hom samples...

      // Using Wakefield's formula for V and ABF
      v = ( casen + controln )/( casen*controln*( 2*maf*revmaf-pow(2*maf*revmaf,2) ) );
      r = w/(v + w);
      temp_abf = exp( pow(temp_Z,2)*r/2 )*sqrt(1-r);
      if( isinf( temp_abf ) ) temp_abf = max_abf;
      (pair_pvec[u])->abf = temp_abf;
      }
    else (pair_pvec[u])->abf = 0;

    cumul_abf += (pair_pvec[u])->abf;
    (pair_pvec[u])->credset = "TRUE";
    }
  }

// Need to sort the pairs by ABF
sort( pair_pvec.begin(), pair_pvec.end(), funct_sort_abfp );

// Then calculate the posterior probs and CredSet true/false
bool credset_bool=0;
if( cumul_abf > 0 )
  {
  for(unsigned u=0; u < pvec_size; u++)
    if( (pair_pvec[u])->candidate_Z() )
      {
      (pair_pvec[u])->pprob = (pair_pvec[u])->abf / cumul_abf;
      cumul_pprob += (pair_pvec[u])->pprob;
      if( credset_bool == 1 ) (pair_pvec[u])->credset = "FALSE";
      if( cumul_pprob > credset_limit ) credset_bool = 1;  // The first over thold is still TRUE
      }
  }
}

void funct_abf_wrapper( vector<class_snp_pair>& pair_vec,
  map<long int, double>& maf_map, map< string, double >& rep_zscore_map )
{  // ABFs are calculated after having the full set of pairs, due to TAD boundary backtracking
unsigned vecsize=0;
vector<class_snp_pair*> temp_pair_vec;
vector<class_snp_pair>::iterator pair_it;
set<string> uniq_rep_set;
set<string>::iterator set_it;
multimap< string, class_snp_pair* > pair_mmapp;
multimap< string, class_snp_pair* >::iterator mmapp_it;
pair< multimap< string, class_snp_pair* >::iterator,
  multimap< string, class_snp_pair* >::iterator > mmapp_pair_it;

// We first remove duplicates (exact same cand/rep rsids)
// arising from overlapping expanded TAD boundaries
cout << "Removing duplicate candidate/reporter pairs: we had " << pair_vec.size() << " pairs,\n";
sort( pair_vec.begin(), pair_vec.end(), funct_sort_rsids );  // Sort by pair RSids
pair_it = unique( pair_vec.begin(), pair_vec.end(), funct_unique_rsids );  // Remove duplicates
pair_vec.resize( distance( pair_vec.begin(), pair_it ) );
vecsize = pair_vec.size();
cout << " and after the removal, we have " << vecsize
  << "\nNow calculating the ABFs and CredSets" << endl;

// Create a multimap with RepRsids -> SNP_pairs
for( unsigned u=0; u<vecsize; u++ )
  {  // This vector does not contain self_pairs for ABF calt
  uniq_rep_set.insert( pair_vec[u].snpp->rsname );
  pair_mmapp.insert( pair<string, class_snp_pair*>( pair_vec[u].snpp->rsname, &(pair_vec[u]) ) );
  }

cout << "We have " << uniq_rep_set.size() << " unique reporters, and a map of size of "
  << pair_mmapp.size() << endl;

for( set_it = uniq_rep_set.begin(); set_it != uniq_rep_set.end(); ++set_it )
  {  // Iterate ALL reporters, and calculate their ABFs
  temp_pair_vec.clear();
  mmapp_pair_it = pair_mmapp.equal_range( *set_it );

  for( mmapp_it = mmapp_pair_it.first; mmapp_it != mmapp_pair_it.second; ++mmapp_it )
    temp_pair_vec.push_back( (*mmapp_it).second );

  // ABF/credset calc is passed a pointer vector, to modify the original elements
  funct_calc_abf( temp_pair_vec, maf_map, rep_zscore_map );
  }
}

void funct_sort_and_writeout( vector<class_snp_pair>& pair_vec,
  class_fstream& out, class_fstream& log_out )
{  // Writes out the genotyped/reporter pairs
unsigned write_counter=0;
stringstream ss;
vector<class_snp_pair>::iterator pair_it;

sort( pair_vec.begin(), pair_vec.end() );  // Sort by ascending Z-score pvals

for( pair_it=pair_vec.begin(); pair_it!=pair_vec.end(); ++pair_it )
  if( (*pair_it).candidate_Z() )
    {  // Writeout - only true candidates with candZ > repZ
    (*pair_it).writeout( out );  // Write standard candidate-rep pairs
    write_counter++;
    }
// Problem of this is that we won't write many control pairs?

ss.str("");
ss << "Writen " << write_counter << " causal-reporter pairs" << endl;
cout << ss.str(); log_out.file << ss.str();
}

/*   =====================
 *   === MAIN FUNCTION ===
 *   =====================
 */

int main( int argc, char* argv[] )
{
unsigned chrn=0, nonhw_count=0;
char header[65536];
string inpath="", ginfopath="", genotypath="",
  outpath="", hi_name="", dummy="", rnadummy="",
  chrt="", strheader="", input_chr="";
stringstream ss;

class_snp snp;
class_genotype_phased genotype;

vector<string> rnaseq_vec, genoty_vec;  // Contains the names of all Samples analised, from the header of trimmed
vector< int > genoty2rnaseq_indx_vec;
vector< pair<double,double> > rndzval_distr;
vector<double> rndzvec(2);  // There are only 2 fields in the control distr: Zscore, density
list<class_snp> repsnp_list;  // Is a list since SNP object is BIG
list<class_snp*> tad_snpps;
vector< vector<class_region> > tad_vec(22);
vector<class_region>::iterator tad_it;
map<string,string> names2gsymb_map;
map<string,string>::iterator map_it;
map<string, double> r2_map, rep_zscore_map;
map<long int, double> maf_map;
set<string> nonhw_rsid_set;
vector<class_genotype_phased> phcand_vec, tphcand_vec;
vector<class_snp_pair> rep_vec, pair_vec;

class_fstream data_in, rnd_zval_in, zscore_in, genoty_in,
  r2_in, maf_in, hwe_in, tads_in,
  reporters_out, credset_out, log_out;

srand( unsigned( time( NULL ) ) );  // To randomly assing direction to AR == 0.5, and create the random controls

if(!debug) cout << endl << "Starting ALLELIC IMBA CREDSET" << endl;
else { cout << endl << "Starting ALLELIC IMBA CREDSET -- in DEBUG MODE --" << endl; }
{  // Obre totes les entrades i sortides
size_t p1=0, p2=0, p3=0;
unsigned field=0;

if( debug )
  {
  inpath = "/home/ender/projects/t2dsys/case/zscore_out/";
  genotypath = "/home/ender/projects/t2dsys/genotypes/";  // t2dsystems_chr22.vcf";
  ginfopath = "/home/ender/Dades/Genomic_Info/";
  outpath = "/home/ender/projects/t2dsys/case/candidates_full/";
  input_chr = testchr;
  }
else if( !debug && argc != 6 ) { blob_of_text(); return 0; }
else if( !debug && argc == 6 )
  {
  inpath     = argv[1];
  genotypath = argv[2];
  ginfopath  = argv[3];
  outpath    = argv[4];
  input_chr  = argv[5];
  }

if( !data_in.open    ( inpath + "Allelic_Imba_trimmed.txt", "in", "data" ) &&
    !rnd_zval_in.open( inpath + "Allelic_Imba_rnd_Zscore.txt", "in", "control_distr" ) &&  // Same as before - was different? Search for rnd script

    !r2_in.open ( genotypath + "plink/r2_" +input_chr+ ".ld", "in", "r2" ) &&
    !maf_in.open( genotypath + "plink/" +input_chr+ ".frq", "in", "MAF" ) &&
    !hwe_in.open( genotypath + "plink/" +input_chr+ ".hwe", "in", "HWE" ) &&
    !tads_in.open ( ginfopath  + "Joan_domains/TAD-like/iTADs_merging_INS_KCN.bed", "in", "TADs" ) &&  // INS/KCN iTAD (artificially?) broken

    !reporters_out.open ( outpath + "Allelic_Imba_reporters_"+input_chr+".txt", "out" ) &&
    !credset_out.open( outpath + "Allelic_Imba_CSet_"+input_chr+".txt", "out" ) &&
    !log_out.open    ( outpath + "Allelic_Imba_CSet_"+input_chr+".log", "out" ) )
  {
  // Get the number and name of samples from the header of Allelic_Imba_trimmed.txt
  data_in.file.getline(header, 65536); strheader = string(header);
  while( ( p2 = strheader.find( "\t", p1 ) ) != string::npos )
    { field++; if( field % 6 == 0 ) { p3 = strheader.find( "_refc", p1 );
        rnaseq_vec.push_back( strheader.substr( p1, p3-p1 ) ); } p1=p2+1; }

  // *** GLOBAL DEFINITION OF NUM_RNASEQ_SAMPLES!! ***
  num_rnaseq_samples = rnaseq_vec.size();
  // *** GLOBAL DEFINITION OF NUM_RNASEQ_SAMPLES!! ***

  // Remove headers
  zscore_in.file.getline(header, 65536);
  r2_in.file.getline(header, 65536);
  maf_in.file.getline(header, 65536);
  hwe_in.file.getline(header, 65536);

  reporters_out.file << "#RepRsname\tZscore\tZpval\t"
    << "jointAR\tmedianCovg\t"
    << "chr\tpose\tNumRef|AltImba\t"
    << "HetARs\tHetPvals\tHetSamples\tRef|AltCounts\n";

  credset_out.file << "#CausalRsname\tRepRsname\tr2\tmedianCovg\t"
    << "Zscore\tZpval\tRepZscore\t"
    << "ABF\tPProb\t95CredSet\t"
    << "jointAR\tchr\tRepPose\tCausalPose\tDist\t"
    << "tad_ini\ttad_fin\tistad\t"
    << "NumHetSignif\tNumHomSignif\tHomZscore\tHomZPval\t"
    << "HetARs\tHetPvals\tHetSamples\tRef|AltCounts\t"
    << "HomARs\tHomPvals\tHomSamples\n";
  }
else { return 1; }

if( !genoty_in.open ( genotypath + "t2dsystems_"+input_chr+".vcf", "in", input_chr+"_vcf" ) )
  {
  genoty_in.file.getline(header, 65536);
  while( header[0] == '#' && header[1] == '#' )
    genoty_in.file.getline(header, 65536);  // Skip all the useless vcf comment lines

  // Gets the names of the genotyped samples in the file (and assumes the same order for the rest of the chr files!)
  strheader = string(header); field=0; p1=0;
  while( ( p2 = strheader.find( "\t", p1 ) ) != string::npos )
    { field++; if( field > 9 ) {
        genoty_vec.push_back( strheader.substr( p1, p2-p1 ) ); } p1=p2+1; }
  genoty_vec.push_back( strheader.substr( p1 ) );

  genoty_in.close();
  }
else return 1;

// *** GLOBAL DEFINITION OF GENOTY_SAMPLES!! ***
num_genoty_samples = genoty_vec.size();
// *** GLOBAL DEFINITION OF GENOTY_SAMPLES!! ***

cout << "We have " << num_rnaseq_samples << " RNAseq and "
  << num_genoty_samples << " genotyped samples" << endl;

genoty2rnaseq_indx_vec.resize( num_rnaseq_samples );
bool gjfound;
for( unsigned k=0; k<num_rnaseq_samples; k++ )
  {  // Find genoty-rnaseq matching names
  gjfound=0;
  rnadummy = rnaseq_vec[k].substr( 0, rnaseq_vec[k].find( "low" ) );  // To accept HIXXlow entries

  for( int j=0; j<(int)num_genoty_samples; j++ )
    if( genoty_vec[j] == rnadummy )
      { genoty2rnaseq_indx_vec[ k ] = j; gjfound=1; }
  if( !gjfound )
    { genoty2rnaseq_indx_vec[ k ] = -1;
      cout << ">> Warning! Could not find a matching genotype for RNAseq " << rnaseq_vec[k] << ". <<" << endl; }
  }
}

cout << "Loading r2, MAF, nonHW, TADS, BenjH and Z-control distr" << endl;

funct_load_r2( r2_in, r2_map );  // Loads the r2 pair values
funct_load_maf( maf_in, maf_map );  // Loads the MAF values
funct_load_hwe( hwe_in, nonhw_rsid_set );  // Loads the nonHWE blacklist
funct_load_tads( tads_in, tad_vec );
for(chrn=0; chrn<22; chrn++)  // Sort
  sort( tad_vec[chrn].begin(), tad_vec[chrn].end() );

// Loads the random distribution controls - updated to sort-smooth
funct_read_rndzscore_and_sort( rnd_zval_in, rndzval_distr );

cout << "Loaded: " << r2_map.size() << " r2 pairs, "
  << maf_map.size() << " MAF, " << nonhw_rsid_set.size()
  << " nonHW rsids, and\n "
  << rndzval_distr.size() << " control distr values\n";

cout << "Reading SNP data" << endl;
chrt=""; snp.read( data_in );
while( !data_in.file.eof() )
  {  // Reads all the SNP data
  if( chrt != snp.chr )
    {
    chrt=snp.chr;
    for(chrn=0; chrn<22; chrn++)
      if( cromosomes[chrn] == chrt ) break;
    }

  if( cromosomes[chrn] == input_chr &&
      snp.num_hets >= min_num_hets )  // Even though they are already prefiltered
    {
    if( nonhw_rsid_set.find( snp.rsname ) == nonhw_rsid_set.end() )
      {  // NonHW filter
      snp.calc_median_coverage();  // Calculates the median het coverage

      repsnp_list.push_back(snp);  // Pushback here since the prev function generates permanent values
      }
    else nonhw_count++;
    }

  snp.read( data_in );
  }
repsnp_list.sort();
cout << "Loaded " << repsnp_list.size() << " SNP data, and discarded "
  << nonhw_count << " SNPs due to nonHW eq" << endl;


// Read the genotypes, and for every genotype,
// calculate all chr Zscores for the hets,
// and writeout the significant reporter pairs

// This assumes the genotype file is sorted by chr but without assuming chr ordering
long fpos;
for( chrn=0; chrn<22; chrn++ )
  if( cromosomes[chrn] == input_chr )
    {  // Loads all genotypes in the RAM and processes the Zscores and overlaps
    if( !genoty_in.open ( genotypath + "t2dsystems_" + cromosomes[chrn] + ".vcf", "in", cromosomes[chrn]+"_vcf" ) )
      {  // Open the phased genotype file, initiate chromosome dependent stuff
      cout << "Reading genotypes for chr" << chrn+1 << " and searching all possible Z-score pairs" << endl;
      genoty_in.file.getline(header, 65536);
      while( header[0] == '#')  // Skip all the useless vcf comment lines
        { fpos = genoty_in.file.tellp();
        genoty_in.file.getline(header, 65536); }
      genoty_in.file.seekp(fpos);  // Rewind 1 line

      // Load reporter SNPs of the TAD, initiate pointers and load r2 map for the chromosome
      tad_it = tad_vec[chrn].begin();
      funct_load_tad_snps( (*tad_it), repsnp_list, tad_snpps );

      phcand_vec.clear();
      }
    else return 1;

    genotype.read( genoty_in );  // MODIFIED CLASS READ to parse vcf format!!
    while( !genoty_in.file.eof() )
      {  // We process each genotyped SNP vs all the reporters of that chr
      if( genotype.alt.find(",") == string::npos &&  // Skip SNPs with 2+ alternate alleles
          nonhw_rsid_set.find( genotype.rsname ) == nonhw_rsid_set.end() )
        {
        if( tad_it != tad_vec[chrn].end() && genotype.pose > (*tad_it).fin + tad_leeway )
          {  // At the end of the TAD, look for all causal-reporter pairs

          // Calculation of all Z-scores and significance occurs HERE!!
          if( !phcand_vec.empty() )
            funct_search_pairs_and_zscores( phcand_vec, tad_snpps, tad_it,
              genoty2rnaseq_indx_vec, r2_map, rndzval_distr,
              rnaseq_vec, rep_vec, rep_zscore_map, pair_vec );

          // Clear all phcand_vec elements that do not belong in the next potentially overlapping TAD
          tphcand_vec.clear();
          for( unsigned uu=0; uu<phcand_vec.size(); uu++ )
            if( phcand_vec[uu].pose >= (*tad_it).fin - tad_leeway )
              tphcand_vec.push_back( phcand_vec[uu] );  // Keep the SNPs near the rightmost border of ending TAD
          phcand_vec = tphcand_vec; tphcand_vec.clear();

          ss.str(""); ss << "Processed " << pair_vec.size() << " causal-reporter pairs so far, in chr" << chrn+1;
          cout << ss.str() << endl; log_out.file << ss.str() << "\n";

          // Loads the SNPs that are located within the next TAD, for pair matching
          while( tad_it != tad_vec[chrn].end() && (*tad_it).fin + tad_leeway < genotype.pose ) ++tad_it;  // Next TAD
          if( tad_it != tad_vec[chrn].end() )
            funct_load_tad_snps( (*tad_it), repsnp_list, tad_snpps );  // And loads the reporter SNPs
          else break;  // End of last TAD, close the chr and to the next one
          }

        if( genotype.pose >= (*tad_it).ini - tad_leeway )
          {
          phcand_vec.push_back( genotype );  // Stores the genotype of the causal SNP with all others in the TAD
          }
        }
      genotype.read( genoty_in );
      }
    genoty_in.close();

    if( !phcand_vec.empty() && !tad_snpps.empty() )
      {  // End of chr pass
      if( tad_it == tad_vec[chrn].end() ) --tad_it;
      cout << "End of CHROM, evaluation of " << phcand_vec.size()
        << " phased vec size" << endl;
      funct_search_pairs_and_zscores( phcand_vec, tad_snpps, tad_it,
        genoty2rnaseq_indx_vec, r2_map, rndzval_distr,
        rnaseq_vec, rep_vec, rep_zscore_map, pair_vec );
      }
    }
cout << "Done with all pairs! Loaded " << pair_vec.size()
  << " pairs and " << rep_vec.size() << " reporters to the vectors"
  << endl;

cout << "Reporter writeout" << endl;
funct_write_reporters( rep_vec, reporters_out );
// Writes ALL, w/o caring about significance

cout << "Now calculating the 95% credible sets" << endl;
funct_abf_wrapper( pair_vec, maf_map, rep_zscore_map );
// We only consider candZ > repZ as candidates, the rest are OUT

cout << "Writing out all candidate-reporter pairs" << endl;
funct_sort_and_writeout( pair_vec, credset_out, log_out );


ss.str(""); ss << "All done!" << endl;
cout << ss.str(); log_out.file << ss.str();


data_in.close(); rnd_zval_in.close();
tads_in.close(); hwe_in.close(); r2_in.close(); maf_in.close();
credset_out.close(); reporters_out.close(); log_out.close();

//cin.get();
}





















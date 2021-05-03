/*
 * Allelic Imba: Filter and Unify
 *
 * Reads all dbSNP142_quantification.txt files, and a .vcf genotype file,
 * with the genotypes in columns in any order (automatic matching w given RNAseq sample names),
 * and filters the data and creates a single file with all the imba information.
 *
 * Last modified: 08/12/16
 */

#include "../includes.h"
#include "../class_fstream.h"
#include "../allelic_imba_common.h"  // Contains the functions to calculate imba, covg criteria

using namespace std;

bool debug = 1;  // DEBUG FLAG
unsigned cloncounter=0;  // Counter of discarded events due to clonality
const unsigned min_gtype_consist = 3;  // Min number of samples to check genotype consistency with RNAseq
const double benj_fdr = 0.01;  // 1% FDR Benjamini-Hochberg threshold per sample

/*
 *   CLASSES
 */

double binom_pval( const pair<int, int>& pr, const double p=0.5 );
double get_chisq_pval( vector<double>& bpval_vec );

class class_bed
{  // Bed3 of the ENCODE forbidden regions
public:
  string chr;
  long int ini, fin;

  void read( class_fstream& in )
    { in.file >> chr >> ini >> fin; }
};

class class_sampleSnp
{  // Conte la info del SNP en aquesta mostra, la quantificacio d'allelic imbalance i el binom pval
public:
  string chr, ref, alt, err, rsname, genoty;
  long int pose, refcount, altcount, errcount;
  double raw_imba, pval, delta;

  void read_and_calculate( class_fstream& in )
    {  // Reads the read counts and calculates the allelic ratio and binom p-value for this sample
    string dummy;
    in.file >> chr >> pose >> ref >> refcount >> alt >> altcount >> err >> errcount
      >> dummy >> rsname;
    if( refcount+altcount > 0 )
      raw_imba = refcount/(double)(refcount+altcount);
    else raw_imba = -1;

    pval = -1; delta = -100;
    genoty = "";  // So that when the values are added to the sampleSnp_vec they dont get removed
    }

  void clean()
    { chr = ""; ref = ""; alt = ""; rsname = ""; pose = 0; }  // To free RAM from unnecessary values

  void zero_all()
    { refcount = -1; altcount = -1; raw_imba = -1; pval = -1; }  // To kill this sampleSNP

  void clonality_delta( class_sampleSnp& clonSnp )
    {  // Checks if there are enough non-clonal reads and calculates the delta
    if( clonSnp.refcount + clonSnp.altcount < min_num_clon_reads )
      { zero_all(); cloncounter++; }
    else { delta = raw_imba - clonSnp.refcount/(double)( clonSnp.refcount + clonSnp.altcount ); }
    }

  void calc_pval( double expected_ratio )
    {  // Calculates the binomial pvalue with the given expected allelic ratio
    if( genoty.substr(0,3) == "Het" )  // Works for Het/Het1/Het2
      pval = binom_pval( pair<int,int>( refcount, altcount ), expected_ratio );
    else pval = -1;
    }

  void correct_AR( double expected_ratio )
    {  // Calculates the corrected ARs so that they are centered at 0.50 .
       // "Raw" ARs can still be obtained by dividing the raw reads
    double imba;
    if( refcount+altcount > 0 ) imba = refcount/(double)(refcount+altcount);
    else imba = -1;

    if( imba > -1 )
      { if( imba <= expected_ratio ) raw_imba = imba*0.50/expected_ratio;
        else raw_imba = 1-( (1-imba)*0.50/(1-expected_ratio) ); }
    }

  bool good_covg()  // Checks if good coverage
    { return refcount + altcount >= min_total_count; }

  bool good_het()
    {  // If this sample is Het, with good coverage and AR
    return genoty.substr(0,3) == "Het" && good_covg() &&
      raw_imba > imba_het_threshold && raw_imba < 1-imba_het_threshold;
    }

  bool good_ref()
    { return genoty == "Ref" && good_covg() && raw_imba >= 1-imba_het_threshold; }
  bool good_alt()
    { return genoty == "Alt" && good_covg() && raw_imba > -1 && raw_imba <= imba_het_threshold; }
};

class class_compiledSnp
{  // Conte la info de totes les mostres per aquest SNP
public:
  unsigned num_good_hets;
  bool good_genotype;
  string chr, ref, alt, rsname, snp_pair;
  long int pose;
  vector<class_sampleSnp> sampleSnp_vec;  // Vector of RNAseq'd samples SNP info

  void init( class_sampleSnp& sp )
    {  // Initialization
    class class_sampleSnp s;

    s.refcount=-1; s.altcount=-1; s.errcount=-1; s.raw_imba=-1; s.pval=-1; // s.genoty="";
    sampleSnp_vec.clear(); sampleSnp_vec.resize( num_rnaseq_samples, s );  // Initialized to empty

    chr = sp.chr; ref = sp.ref; alt = sp.alt; rsname = sp.rsname; pose = sp.pose; snp_pair = ref+alt;
    }

  bool usable_snp()
    {  // Checks if it has 3+ Het SNPs with min_num_reads+, reasonable ARs, and consistent genotypes
    unsigned num_ngt=0;
    double refm, hetm, altm;
    vector<double> refvec, hetvec, altvec;  // Contain ALL Ref, Het and Alt samples AR

    num_good_hets=0; good_genotype=0;

    for( unsigned k=0; k<num_rnaseq_samples; k++ )
      {
      if( sampleSnp_vec[k].genoty == "Ref" && sampleSnp_vec[k].raw_imba > -1 )
        refvec.push_back( sampleSnp_vec[k].raw_imba );
      else if( sampleSnp_vec[k].genoty.substr(0,3) == "Het" && sampleSnp_vec[k].raw_imba > -1 )
        {
        if( sampleSnp_vec[k].good_het() ) num_good_hets++;
        hetvec.push_back( sampleSnp_vec[k].raw_imba );
        }
      else if( sampleSnp_vec[k].genoty == "Alt" && sampleSnp_vec[k].raw_imba > -1 )
        altvec.push_back( sampleSnp_vec[k].raw_imba );
      else if( sampleSnp_vec[k].genoty == "" ) num_ngt++;
      }

    // Genotype consistency check by checking the median imbas versus reasonable ARs
    if( refvec.size() >= min_gtype_consist ) refm = funct_get_median( refvec ); else refm = -1;
    if( hetvec.size() >= min_gtype_consist ) hetm = funct_get_median( hetvec ); else hetm = -1;
    if( altvec.size() >= min_gtype_consist ) altm = funct_get_median( altvec ); else altm = -1;

    if( num_ngt != num_rnaseq_samples &&
        ( hetm > imba_het_threshold && hetm < 1-imba_het_threshold ) &&
        ( refm > 1-imba_het_threshold || refm < 0 ) &&  // To allow for non Ref samples
        ( altm < imba_het_threshold ) )
      good_genotype = 1;

    if( num_good_hets >= min_num_hets && good_genotype == 1 ) return 1;
    return 0;
    }

  void writeout( class_fstream& full_out )
    {  // Writes out (d'oh)
    string gty;

    full_out.file << chr << "\t" << pose << "\t" << rsname << "\t" << ref << "\t" << alt;
    for( unsigned k=0; k<num_rnaseq_samples; k++ )
      {  // Prints everything
      if( sampleSnp_vec[k].genoty == "" ) gty = "NGT";
      else gty = sampleSnp_vec[k].genoty;
      full_out.file << "\t" << sampleSnp_vec[k].refcount << "\t" << sampleSnp_vec[k].altcount << "\t"
        << sampleSnp_vec[k].errcount << "\t" << sampleSnp_vec[k].raw_imba << "\t"
        << sampleSnp_vec[k].pval << "\t" << gty;
      }
    full_out.file << "\n";
    }
};

class class_count_stats
{  // Contains the ref/alt counts and vector of ARs
public:
  long int refc, altc;
  double raw_AR, mean_AR, median_AR;
  vector< double > AR_vec;

  class_count_stats() : raw_AR(-1), mean_AR(-1), median_AR(-1) {}
};

/*
 *   FUNCIONS SECUNDARIES
 */

bool funct_sort_bed( const class_bed &b1, const class_bed &b2 )
{ return b1.ini < b2.ini; }

void blob_of_text()
{  // A blob of text for when the wrong arguments are passed
cout << "\nThis program reads a list of SNP read count quantification files, "
  << "and filters and unifies them.\n\n"
  << "It requires 5 arguments:\n"
  << " 1) The general path to the folders where the files are located, ending in /. "
  << "Ie /project/jferrer/imperial/Fastq_Data/\n"
  << " 2) A comma-separated list of folder names where the dbSNP142_quantification.txt are located. "
  << "Ie HI12,HI14oxf,HI15,HI17 which will use the following files:\n"
  << "/project/jferrer/imperial/Fastq_Data/HI12/merged_out/dbSNP142_quantification.txt\n"
  << "/project/jferrer/imperial/Fastq_Data/HI12/merged_nonclonal/dbSNP142_quantification.txt\n"
  << "/project/jferrer/imperial/Fastq_Data/HI14oxf/merged_out/dbSNP142_quantification.txt\n"
  << "/project/jferrer/imperial/Fastq_Data/HI14oxf/merged_nonclonal/dbSNP142_quantification.txt , etc.\n"
  << "The program assumes that the files are named dbSNP142_quantification.txt and the folders "
  << "merged_out and merged_nonclonal.\n"
  << " 3) The path to a file containing the sample genotypes, which must be in .vcf format "
  << "(phased or unphased), and which must contain the same names in the header as the ones "
  << "provided in 2).\n"
  << " 4) The path to a bed3 file containing regions to be excluded from the analysis. "
  << "The first line of the file will always be skipped, so it should be a header. Ie\n"
  << "# From  http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability\n"
  << "chr1 564449 570371\nchr1 724136 727043\n"
  << " 5) The path where the output files will be located, ending in /. "
  << "Ie /project/jferrer/imperial/Fastq_Data/Results/\n"
  << "\nFor example, a valid call would be\n"
  << "$ ./allelic_imba_filter_and_unify /project/jferrer/imperial/Fastq_Data/ HI12,HI14oxf,HI15,HI17 "
  << "/project/jferrer/imperial/Fastq_Data/Genotypes/Islet_SNP_genotypes.vcf "
  << "/project/jferrer/imperial/Fastq_Data/Genotypes/ENCODE_blacklisted.bed "
  << "/project/jferrer/imperial/Fastq_Data/Results/\n"
  << "\nIf any of the above expectations are not met, the program will most probably "
  << "simply crash without much explanation or get stuck in an infinite loop.\n"
  << "\nThere are four output files:\n"
  << " - Allelic_Imba_compiled.txt , which contains, for all SNPs in all samples, "
  << "the reference and alternate read counts, the raw allelic ratio and their genotype.\n"
  << " - Allelic_Imba_compiled_nonclonal.txt , which contains the same information but "
  << "after removal of clonal reads.\n"
  << " - Allelic_Imba_trimmed.txt , which contains, for only those SNPs with usable Het values in 3+ samples, "
  << "the reference and alternate read counts, the corrected allelic ratio and p-value and their genotype.\n"
  << " - Allelic_Imba_AR_biases.txt , which contains the mean and median allelic ratio biases "
  << "per sample per dinucleotide.\n\n"
  << "This software is provided without any warranty whatsoever, so use it at your own risk.\n";
}

double binom_prob( const int m, const int n, const double p )
{  // Given M good and N fail with probs p and q, returns the binom probability
double q = 1-p;

double temp = lgamma(m + n + 1.0);
temp -= lgamma(n + 1.0) + lgamma(m + 1.0);
temp += m*log(p) + n*log(q);

return exp(temp);
}

double binom_pval( const pair<int, int>& pr, const double p )  // p=0.5 by default, already declared
{  // The 2-tailed pval is the cumulative probability of the most extreme cases
int m=pr.first, n=pr.second, i, j;
double epsilon = 1e-7, b, limit = binom_prob( m, n, p ), cumul=0;

if( m+n<1 ) return -1;  // Too few reads

i=0, j=m+n;
while( ( b = binom_prob( i, j, p ) ) <= limit && i <= m+n )
  { cumul += b; i++; j--; }

if( cumul < 1-epsilon ) {  // If the given point is the peak, 1st pass will ==1. Otherwise, this executes
  i=m+n, j=0;
  while( ( b = binom_prob( i, j, p ) ) <= limit && i >= 0 )
    { cumul += b; i--; j++; } }

if( cumul > 1+epsilon ) { cout << "Binom pval error: pval of " << cumul << " with " << m << ", " << n
  << " and p " << p << "!\n"; cin.get(); }

return cumul;
}

void funct_genoty_the_snps( class_genotype& genoty, class_compiledSnp& compSnp,
  vector<int>& genoty2rnaseq_indx_vec )
{  // Genotypes the sampleSnps of compiledSnp
map< string, unsigned >::iterator map_it;

// Using the genoty/rnaseq match index vector
for( unsigned k=0; k<num_rnaseq_samples; k++ )
  {
  if( genoty2rnaseq_indx_vec[k] != -1 )
    compSnp.sampleSnp_vec[k].genoty = genoty.genoty_vec[ genoty2rnaseq_indx_vec[k] ];
  else compSnp.sampleSnp_vec[k].genoty = "NAN";
  }
}

void funct_calc_AR_medians_and_writeout(
  vector< map< string, class_compiledSnp > >& compile_map,
  vector< string >& rnaseq_vec, vector< string >& snp_pairs,
  vector< vector< class_count_stats > >& stats_vec, class_fstream& out )
{  // Calculates the per sample per dint allelic ratio biases and outputs them
bool warn;
unsigned vecsize;
long int refc, altc;
double m;
vector< unsigned > good( num_rnaseq_samples,0 ), total( num_rnaseq_samples,0 );
map< string, class_compiledSnp >::iterator comp_it;

// Adds up all usable het ref and alt reads per sample and/or snp
for(unsigned chrn=0; chrn<22; chrn++)
  for( comp_it = compile_map[chrn].begin(); comp_it != compile_map[chrn].end(); ++comp_it )
    if( (*comp_it).second.usable_snp() )  // Contains 3+ good hets
      for( unsigned snp_indx=0; snp_indx<12; snp_indx++ )  // 12 dinucleotide permutations
        if( (*comp_it).second.snp_pair == snp_pairs[snp_indx] )  // Separates SNPs by allelic pairs
          for( unsigned sample=0; sample<num_rnaseq_samples; sample++ )
            {
            if( (*comp_it).second.sampleSnp_vec[ sample ].good_het() )  // and Het is usable
              {
              refc = (*comp_it).second.sampleSnp_vec[ sample ].refcount;
              altc = (*comp_it).second.sampleSnp_vec[ sample ].altcount;

              // Compiles the number of ref and alt reads per sample per snp
              stats_vec[ sample ][ snp_indx ].refc += refc;
              stats_vec[ sample ][ snp_indx ].altc += altc;
              stats_vec[ sample ][ snp_indx ].
                AR_vec.push_back( refc/(double)(refc+altc) );

              good[sample]++; total[sample]++;
              }
            else if( (*comp_it).second.sampleSnp_vec[ sample ].good_ref() ||
              (*comp_it).second.sampleSnp_vec[ sample ].good_alt() ) { good[sample]++; total[sample]++; }
            else if ( (*comp_it).second.sampleSnp_vec[ sample ].good_covg() &&
              (*comp_it).second.sampleSnp_vec[ sample ].genoty != "NAN" ) total[sample]++;
            }

// Sanity check comparing RNAseq and genotypes
warn = 0;
for( unsigned sample=0; sample<num_rnaseq_samples; sample++ )
  if( good[sample]/(double)total[sample] < 0.90 ) { warn = 1; break; }
if( warn )
  {
  cout << "\nWARNING: Big discrepancy between RNA and genotype detected - "
    << "please make sure that the genotype file follows the same order than the "
    << "comma separated list of samples. The expected concordance values are >90%, "
    << "values around 50-60% indicate an unpaired sample. The values were:\n";
  for( unsigned sample=0; sample<num_rnaseq_samples; sample++ )
    cout << rnaseq_vec[sample] << "\t" << good[sample]*100/(double)total[sample] << "%\n";
  cout << "\n";
  }

// Per Sample per SNP calculate rawm mean and median ARs and printout
out.file << "#Name\tRawAR\tMeanAR\tMedianAR\n";
for( unsigned sample=0; sample<num_rnaseq_samples; sample++ )
  for( unsigned snp_it=0; snp_it<12; snp_it++ )
    {
    stats_vec[sample][snp_it].raw_AR =
      stats_vec[sample][snp_it].refc/(double)( stats_vec[sample][snp_it].refc +
        stats_vec[sample][snp_it].altc );

    m = 0; vecsize = stats_vec[sample][snp_it].AR_vec.size();
    for( unsigned k=0; k<vecsize; k++ ) m += stats_vec[sample][snp_it].AR_vec[k];
    stats_vec[sample][snp_it].mean_AR = m/(double)vecsize;

    stats_vec[sample][snp_it].median_AR = funct_get_median( stats_vec[sample][snp_it].AR_vec );

    out.file << rnaseq_vec[sample] << "_" << snp_pairs[snp_it] << "\t"
      << stats_vec[sample][snp_it].raw_AR << "\t"
      << stats_vec[sample][snp_it].mean_AR << "\t" << stats_vec[sample][snp_it].median_AR << "\n";
    }
}

/*   =====================
 *   === MAIN FUNCTION ===
 *   =====================
 */

int main( int argc, char* argv[] )
{
unsigned chrn=0, vecsize;
char header[65536];
string path, hilist, genotypath, blackpath, outpath, rnadummy,
  standard_path, nonclonal_path, chrt="", cromosomes[22] =
  {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
   "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
   "chr20","chr21","chr22"},
  nucleotides[4]={"A","C","G","T"};
stringstream ss;

class_bed bed;
class_sampleSnp sampleSnp, clonSnp;
class_genotype genotype;
class_compiledSnp compiledSnp;

vector< int > genoty2rnaseq_indx_vec;
vector< double > delta_vec;
vector< vector< double > > benj_vec;
vector< string > rnaseq_vec, genoty_vec, snp_pairs;
vector< vector<class_bed> > forbidden_vec(22);
vector<class_bed>::iterator bed_it;
vector< map< string, class_compiledSnp > > compile_map(22), nonclon_map(22);
map< string, class_compiledSnp >::iterator comp_it, noncl_it;
vector< vector< class_count_stats > > AR_bias_supervec;

class_fstream* hi_inp;
class_fstream forbidden_in, genoty_in,
  trim_out, fullinfo_out, fullnoncl_out, AR_bias_out, benj_out;
vector<class_fstream*> snp_invec, nc_invec;

for(unsigned i=0; i<4; i++) for(unsigned j=0; j<4; j++) if(i != j)
  snp_pairs.push_back( nucleotides[i]+nucleotides[j] );  // Creates the 12 dinucleotide pairs

if(!debug) cout << "Starting ALLELIC IMBA FILTER AND UNIFY" << endl;
else cout << "Starting ALLELIC IMBA FILTER AND UNIFY -- in DEBUG MODE --" << endl;
{  // Obre totes les entrades i sortides
if( debug )
  {
  path = "/home/bscuser/project_data/t2dsystems/case/pipeline_tests/";
// firstbatch
//  hilist = "HI12,HI14oxf,HI15,HI16,HI17,HI19,HI20,HI22,HI24,HI25_51,HI26,HI27,HI27oxf,HI28oxf,HI29,HI30,HI31,HI32_51,HI33,HI34,HI37oxf,HI45,HI46,HI57,HI76,HI77,HI81,HI82";
// secondbatch
//  hilist = "HI100,HI101,HI102,HI103,HI111,HI112,HI113,HI115,HI116,HI117,HI118,HI121,HI122,HI123,HI124,HI125,HI126,HI127,HI128,HI79,HI84,HI86,HI87,HI88,HI90,HI91,HI95,HI97,HI99";
// first+thirdbatch
//  hilist = "HI06,HI07,HI12,HI14,HI14oxf,HI15,HI16,HI17,HI18,HI19,HI20,HI22,HI24,HI25_51,HI26,HI27,HI27oxf,HI28oxf,HI29,HI30,HI31,HI32_51,HI33,HI34,HI37oxf,HI44,HI45,HI46,HI57,HI73,HI76,HI77,HI78,HI81,HI82,HI23,HI35,HI37";
// second+thirdbatch
//  hilist = "HI09,HI11,HI13,HI21,HI28,HI79,HI80,HI83,HI84,HI86,HI87,HI88,HI89,HI90,HI91,HI95,HI97,HI99,HI100,HI101,HI102,HI103,HI111,HI112,HI113,HI114,HI115,HI116,HI117,HI118,HI119,HI121,HI122,HI123,HI124,HI125,HI126,HI127,HI128,HI129,HI130,HI131,HI133,HI134,HI135,HI137,HI132,HI138,HI139";
// 1+2+3
//  hilist = "HI06,HI07,HI100,HI101,HI102,HI103,HI111,HI112,HI113,HI114,HI115,HI116,HI117,HI118,HI119,HI12,HI123,HI124,HI125,HI126,HI127,HI128,HI129,HI13,HI130,HI131,HI133,HI134,HI135,HI137,HI14,HI14oxf,HI15,HI17,HI18,HI19,HI20,HI21,HI22,HI24,HI25_51,HI26,HI27,HI28,HI28oxf,HI29,HI30,HI31,HI32_51,HI33,HI34,HI37oxf,HI44,HI45,HI46,HI57,HI73,HI76,HI77,HI78,HI79,HI80,HI81,HI82,HI83,HI84,HI86,HI87,HI88,HI89,HI90,HI91,HI95,HI99";

// ALL samples 1+2+3+4+LG
//hilist = "HI06,HI07,HI100,HI101,HI102,HI103,HI111,HI112,HI113,HI114,HI115,HI116,HI117,HI118,HI119,HI12,HI123,HI124,HI125,HI126,HI127,HI128,HI129,HI13,HI130,HI131,HI133,HI134,HI135,HI137,HI138,HI139,HI14,HI141,HI142,HI143,HI144,HI145,HI147,HI148,HI149,HI14oxf,HI15,HI151,HI152,HI17,HI18,HI19,HI20,HI21,HI22,HI23,HI24,HI26,HI27,HI28,HI28oxf,HI29,HI30,HI31,HI33,HI34,HI37oxf,HI44,HI45,HI46,HI57,HI73,HI76,HI77,HI78,HI79,HI80,HI81,HI82,HI83,HI84,HI86,HI87,HI88,HI89,HI90,HI91,HI95,HI99,ID10,ID12,ID14,ID15,ID16,ID17,ID18,ID19,ID20,ID6,ID7,ID8,ID9,T2D_3,T2D_4,T2D_5,T2D_6,T2D_7,T2D_8,T2D_9,LG_HI1,LG_HI10,LG_HI100,LG_HI101,LG_HI102,LG_HI103,LG_HI104,LG_HI105,LG_HI106,LG_HI107,LG_HI108,LG_HI109,LG_HI110,LG_HI12,LG_HI13,LG_HI14,LG_HI15,LG_HI16,LG_HI19,LG_HI2,LG_HI20,LG_HI21,LG_HI22,LG_HI23,LG_HI24,LG_HI25,LG_HI26,LG_HI3,LG_HI30,LG_HI31,LG_HI32,LG_HI33,LG_HI34,LG_HI35,LG_HI36,LG_HI37,LG_HI38,LG_HI39,LG_HI4,LG_HI40,LG_HI42,LG_HI43,LG_HI44,LG_HI45,LG_HI46,LG_HI48,LG_HI49,LG_HI5,LG_HI50,LG_HI51,LG_HI52,LG_HI53,LG_HI54,LG_HI55,LG_HI56,LG_HI58,LG_HI59,LG_HI6,LG_HI60,LG_HI62,LG_HI64,LG_HI65,LG_HI66,LG_HI67,LG_HI7,LG_HI70,LG_HI71,LG_HI72,LG_HI73,LG_HI75,LG_HI76,LG_HI77,LG_HI78,LG_HI79,LG_HI8,LG_HI80,LG_HI81,LG_HI82,LG_HI84,LG_HI85,LG_HI86,LG_HI87,LG_HI88,LG_HI89,LG_HI9,LG_HI90,LG_HI91";
// 1+2+3+4 samples
//hilist = "HI06,HI07,HI100,HI101,HI102,HI103,HI111,HI112,HI113,HI114,HI115,HI116,HI117,HI118,HI119,HI12,HI123,HI124,HI125,HI126,HI127,HI128,HI129,HI13,HI130,HI131,HI133,HI134,HI135,HI137,HI138,HI139,HI14,HI141,HI142,HI143,HI144,HI145,HI147,HI148,HI149,HI14oxf,HI15,HI151,HI152,HI17,HI18,HI19,HI20,HI21,HI22,HI23,HI24,HI26,HI27,HI28,HI28oxf,HI29,HI30,HI31,HI33,HI34,HI37oxf,HI44,HI45,HI46,HI57,HI73,HI76,HI77,HI78,HI79,HI80,HI81,HI82,HI83,HI84,HI86,HI87,HI88,HI89,HI90,HI91,HI95,HI99,ID10,ID12,ID14,ID15,ID16,ID17,ID18,ID19,ID20,ID6,ID7,ID8,ID9,T2D_3,T2D_4,T2D_5,T2D_6,T2D_7,T2D_8,T2D_9";

// Max debug list
//hilist = "HI06,HI07,HI100,HI101,HI102,HI103,HI111,HI112,HI113,HI114,HI115,HI116,HI117,HI118,HI119,HI12,HI123,HI124,HI125,HI126,HI127,HI128,HI129,HI13,HI130,HI131,HI133,HI134,HI135,HI137,HI138,HI139";

// LG samples
//hilits = "LG_HI1,LG_HI10,LG_HI100,LG_HI101,LG_HI102,LG_HI103,LG_HI104,LG_HI105,LG_HI106,LG_HI107,LG_HI108,LG_HI109,LG_HI110,LG_HI12,LG_HI13,LG_HI14,LG_HI15,LG_HI16,LG_HI19,LG_HI2,LG_HI20,LG_HI21,LG_HI22,LG_HI23,LG_HI24,LG_HI25,LG_HI26,LG_HI3,LG_HI30,LG_HI31,LG_HI32,LG_HI33,LG_HI34,LG_HI35,LG_HI36,LG_HI37,LG_HI38,LG_HI39,LG_HI4,LG_HI40,LG_HI42,LG_HI43,LG_HI44,LG_HI45,LG_HI46,LG_HI48,LG_HI49,LG_HI5,LG_HI50,LG_HI51,LG_HI52,LG_HI53,LG_HI54,LG_HI55,LG_HI56,LG_HI58,LG_HI59,LG_HI6,LG_HI60,LG_HI62,LG_HI64,LG_HI65,LG_HI66,LG_HI67,LG_HI7,LG_HI70,LG_HI71,LG_HI72,LG_HI73,LG_HI75,LG_HI76,LG_HI77,LG_HI78,LG_HI79,LG_HI8,LG_HI80,LG_HI81,LG_HI82,LG_HI84,LG_HI85,LG_HI86,LG_HI87,LG_HI88,LG_HI89,LG_HI9,LG_HI90,LG_HI91";

// Minimal hilist
  hilist = "HI100,HI102,HI103";

  genotypath = "/home/bscuser/project_data/t2dsystems/case/pipeline_tests/minimal_genotypes.vcf";
  blackpath  = "/home/bscuser/Genomic_Info/ENCODE_BROAD_blacklisted.bed";
  outpath    = path + "original/";
  }
else if( !debug && argc != 6 ) { blob_of_text(); cout << "Terminating.\n"; return 1; }
else if( !debug && argc == 6 )
  {  // Arguments from command line
  path       = argv[1];
  hilist     = argv[2];
  genotypath = argv[3];
  blackpath  = argv[4];
  outpath    = argv[5];

  cout << "Given paramters were:\n";
  for(unsigned u=1; u<6; u++) cout << u << ") " << argv[u] << "\n";
  cout << "\n";
  }
standard_path  = "merged_out/";  // Non optional at the moment
nonclonal_path = "merged_nonclonal/";

// Reads the comma separated hilist and puts it in a string vector
size_t p1=0,p2=0;
while( ( p2 = hilist.find( ",", p1 ) ) != string::npos )
  { rnaseq_vec.push_back( hilist.substr( p1, p2-p1 ) ); p1 = p2+1; }
rnaseq_vec.push_back( hilist.substr( p1, string::npos ) );  // Last element
if( rnaseq_vec.back() == "" ) rnaseq_vec.pop_back();  // Solves strings ending in comma

// *** GLOBAL DEFINITION OF num_rnaseq_samples HERE ***
num_rnaseq_samples = rnaseq_vec.size();
// *** GLOBAL DEFINITION OF num_rnaseq_samples HERE ***

// And a few vector sizes that depend on it
AR_bias_supervec.resize(
  num_rnaseq_samples, vector< class_count_stats >(12) );  // 12 dinucleotide combinations
benj_vec.resize( num_rnaseq_samples );

for( unsigned u=0; u<num_rnaseq_samples; u++ )
  {  // Pointer vecs to all quantification files, standard and non-clonal
  try { hi_inp = new class_fstream; }  // Using NEW and DELETE so that the file does not go out of scope
  catch (bad_alloc& ba) { cout << "bad_alloc caught: " << ba.what() << "\n"; }

  if( !hi_inp->obre( path + rnaseq_vec[u] + "/" + standard_path + "dbSNP142_quantification.txt",
      "in", rnaseq_vec[u] + "_standard") )
    snp_invec.push_back( hi_inp );
  else return 1;

  // Now non-clonal
  try { hi_inp = new class_fstream; }  // Using NEW and DELETE so that the file does not go out of scope
  catch (bad_alloc& ba) { cout << "bad_alloc caught: " << ba.what() << "\n"; }

  if( !hi_inp->obre( path + rnaseq_vec[u] + "/" + nonclonal_path + "dbSNP142_quantification.txt",
      "in", rnaseq_vec[u] + "_nonclonal" ) )
    nc_invec.push_back( hi_inp );
  else return 1;
  }

if( !genoty_in.obre   ( genotypath, "in", "genotype" ) &&
    !forbidden_in.obre( blackpath,  "in", "blacklist" ) &&

    !AR_bias_out.obre   ( outpath + "Allelic_Imba_AR_biases.txt", "out" ) &&  // Per-sample per-dinucleotide allelic ratio biases, for posterior correction
    !fullinfo_out.obre  ( outpath + "Allelic_Imba_compiled.txt", "out" ) && // All info, pre filtering
    !fullnoncl_out.obre ( outpath + "Allelic_Imba_compiled_nonclonal.txt", "out" ) && // Non-clonal info, pre filtering
    !trim_out.obre      ( outpath + "Allelic_Imba_trimmed.txt", "out" ) &&  // Trimmed results to only 3+ valid heterozgous SNPs
    !benj_out.obre      ( outpath + "Allelic_Imba_Benjamini.txt", "out" ) )  // Per-sample Benjamini Hockberg thresholds
  {
  string strheader;
  size_t p1=0, p2=0;
  unsigned field=0;

  genoty_in.file.getline(header, 65536);
  while( header[0] == '#' && header[1] == '#' ) genoty_in.file.getline(header, 65536);
  strheader = string(header); field=0; p1=0;
  while( ( p2 = strheader.find( "\t", p1 ) ) != string::npos )
    { field++; if( field > 9 ) {
        genoty_vec.push_back( strheader.substr( p1, p2-p1 ) ); } p1=p2+1; }
  genoty_vec.push_back( strheader.substr( p1 ) );

  forbidden_in.file.getline(header, 65536);

  ss << "#chr\tposition\trsname\tref\talt";
  for( unsigned k=0; k<num_rnaseq_samples; k++ )
    ss << "\t" << rnaseq_vec[k] << "_refc\t" << rnaseq_vec[k] << "_altc\t"
      << rnaseq_vec[k] << "_errc\t" << rnaseq_vec[k] << "_imba\t"
      << rnaseq_vec[k] << "_pval\t" << rnaseq_vec[k] << "_genoty";
  ss << "\n";
  fullinfo_out.file << ss.str(); trim_out.file << ss.str();
  }
else return 1;

// *** GLOBAL DEFINITION OF NUM_GENOTY_SAMPLES!! ***
num_genoty_samples = genoty_vec.size();
// *** GLOBAL DEFINITION OF NUM_GENOTY_SAMPLES!! ***

// Now calculate the RNAseq/genoty index concordance vector
cout << "We have " << num_rnaseq_samples << " RNAseq and " << num_genoty_samples
  << " genotype samples\n";

genoty2rnaseq_indx_vec.resize( num_rnaseq_samples );
bool gjfound;
for( unsigned k=0; k<num_rnaseq_samples; k++ )
  {
  gjfound=0;
  rnadummy = rnaseq_vec[k].substr( 0, rnaseq_vec[k].find( "low" ) );  // To accept HIXXlow entries

  for( int j=0; j<(int)num_genoty_samples; j++ )
    if( genoty_vec[j] == rnadummy )
      { genoty2rnaseq_indx_vec[ k ] = j; gjfound=1; }
  if( !gjfound )
    { genoty2rnaseq_indx_vec[ k ] = -1;
      cout << ">> Warning! Could not find a matching genotype for " << rnaseq_vec[k] << ". <<\n"; }
  }
}

cout << "Reading the forbidden regions\n";
chrt=""; bed.read( forbidden_in );
while( !forbidden_in.file.eof() )
  {  // Loads the forbidden regions
  if( bed.chr != chrt )
    {
    chrt = bed.chr;
    for(chrn=0; chrn<22; chrn++) if( chrt == cromosomes[chrn] ) break;
    }

  if(chrn<22) forbidden_vec[chrn].push_back(bed);

  bed.read( forbidden_in );
  }
for(chrn=0; chrn<22; chrn++)  // Sorts the forbidden regions
  sort( forbidden_vec[chrn].begin(), forbidden_vec[chrn].end(), funct_sort_bed );

cout << "Compiling all SNP values over all samples\n";
for(unsigned sample=0; sample<num_rnaseq_samples; sample++)
  {  // Reads all sampleSnps and calculates their imba and pvalue values
  cout << "Working on sample " << rnaseq_vec[sample] << "\n";
  sampleSnp.read_and_calculate( *snp_invec[sample] );
  while( !snp_invec[sample]->file.eof() )
    {
    if( chrt != sampleSnp.chr )
      {  // Chr change
      chrt = sampleSnp.chr;
      for(chrn=0; chrn<22; chrn++) if( cromosomes[chrn] == chrt) break;
      bed_it = forbidden_vec[chrn].begin();
      }

    if(chrn<22)
      {  // Reads the SNPs and puts them and their imba values in the maps
      comp_it = compile_map[chrn].find( sampleSnp.rsname );

      while( bed_it != forbidden_vec[chrn].end() &&
             (*bed_it).fin < sampleSnp.pose ) ++bed_it;

      if( sampleSnp.pose > (*bed_it).fin || sampleSnp.pose < (*bed_it).ini )
        {  // If SNP does NOT overlap a forbidden region
        if( comp_it == compile_map[chrn].end() )
          {  // First time we see the SNP
          compiledSnp.init( sampleSnp );  // Initialization
          sampleSnp.clean();  // Strip repeated values for RAM efficiency
          compiledSnp.sampleSnp_vec[ sample ] = sampleSnp;
          compile_map[chrn][ compiledSnp.rsname ] = compiledSnp;
          }
        else
          {
          sampleSnp.clean();  // Strip repeated values for RAM efficiency
          (*comp_it).second.sampleSnp_vec[ sample ] = sampleSnp;
          }
        }
      }

    sampleSnp.read_and_calculate( *snp_invec[sample] );
    }

//    for(unsigned u=0; u<22; u++)
//      cout << "chr" << u << " has " << compile_map[u].size() << " elements\n"; cin.get();

  }

while( !snp_invec.empty() )
  {  // Removes all new pointers without leaking
  hi_inp = snp_invec.back();
  hi_inp->tanca();  // Closes the associated file
  snp_invec.pop_back();
  delete hi_inp;
  }

// Crosses the genotype info with the sampleSnp info
cout << "Loading the genotypes\n";
chrt=""; genotype.read_vcf( genoty_in );
while( !genoty_in.file.eof() )
  {  // Puja els genotips en un rsname map
  if( chrt != genotype.chr )
    { chrt = genotype.chr; for(chrn=0; chrn<22; chrn++) if( cromosomes[chrn] == chrt) break;
      cout << "Loading chr" << chrn+1 << "\n"; }

  if(chrn<22)
    {
    comp_it = compile_map[chrn].find( genotype.rsname );
    if( comp_it != compile_map[chrn].end() )
      funct_genoty_the_snps( genotype, (*comp_it).second, genoty2rnaseq_indx_vec );
    }

  genotype.read_vcf( genoty_in );
  }

cout << "Writing the compiled noncorrected SNP data\n";
for(chrn=0; chrn<22; chrn++)
  for( comp_it = compile_map[chrn].begin(); comp_it != compile_map[chrn].end(); ++comp_it )
    (*comp_it).second.writeout( fullinfo_out );

// Removes SNPs with less than 15 clonal reads
cout << "Reading all nonclonal values and compiling them\n";
for(unsigned sample=0; sample<num_rnaseq_samples; sample++)
  {  // Reads all sampleSnps and calculates their imba and pvalue values
  cout << "Now working on sample " << rnaseq_vec[sample] << "\n";
  clonSnp.read_and_calculate( *nc_invec[sample] );
  while( !nc_invec[sample]->file.eof() )
    {
    if( chrt != clonSnp.chr )
      {  // Chr change
      chrt = clonSnp.chr;
      for(chrn=0; chrn<22; chrn++) if( cromosomes[chrn] == chrt) break;
      bed_it = forbidden_vec[chrn].begin();
      }

    if(chrn<22)
      {  // Reads the non-clonal data, compiles it in a nc_map and removes the clonal SNPs from compile_map
      while( bed_it != forbidden_vec[chrn].end() &&
             (*bed_it).fin < clonSnp.pose ) ++bed_it;

      comp_it = compile_map[chrn].find( clonSnp.rsname );
      if( ( clonSnp.pose > (*bed_it).fin || clonSnp.pose < (*bed_it).ini ) &&
          ( comp_it != compile_map[chrn].end() ) )  // If SNP does not overlap a forbidden region
        { // Calculates the delta and zeroes the clonSnp if clonal
        (*comp_it).second.sampleSnp_vec[sample].clonality_delta( clonSnp );

        noncl_it = nonclon_map[chrn].find( clonSnp.rsname );
        if( noncl_it == nonclon_map[chrn].end() )
          {  // First time we see the SNP
          compiledSnp.init( clonSnp );  // Initialization
          clonSnp.clean();  // Strip repeated values for RAM efficiency
          compiledSnp.sampleSnp_vec[ sample ] = clonSnp;
          compiledSnp.sampleSnp_vec[ sample ].genoty =
            (*comp_it).second.sampleSnp_vec[ sample ].genoty;  // Adds the genotype info
          nonclon_map[chrn][ compiledSnp.rsname ] = compiledSnp;
          }
        else
          {
          clonSnp.clean();  // Strip repeated values for RAM efficiency
          (*noncl_it).second.sampleSnp_vec[ sample ] = clonSnp;
          (*noncl_it).second.sampleSnp_vec[ sample ].genoty =
            (*comp_it).second.sampleSnp_vec[ sample ].genoty;  // Adds the genotype info
          }
        }
      }

    clonSnp.read_and_calculate( *nc_invec[sample] );
    }
  }

cout << cloncounter << " events have been discarded due to clonality\n"
  << "Now writing the compiled noncorrected non-clonal SNP data\n";
for(chrn=0; chrn<22; chrn++)
  { for( noncl_it = nonclon_map[chrn].begin(); noncl_it != nonclon_map[chrn].end(); ++noncl_it )
      (*noncl_it).second.writeout( fullnoncl_out );
    nonclon_map[chrn].clear(); }
nonclon_map.clear();  // And removes the nonclonal map to free the ram

while( !nc_invec.empty() )
  {  // Removes all new pointers without leaking
  hi_inp = nc_invec.back();
  hi_inp->tanca();  // Closes the associated file
  nc_invec.pop_back();
  delete hi_inp;
  }

// Now calculate the AR deviations and correct them
cout << "Calculating and writing the raw AR biases\n";
funct_calc_AR_medians_and_writeout( compile_map, rnaseq_vec, snp_pairs, AR_bias_supervec, AR_bias_out );

cout << "Correcting the bias and writing the trimmed SNP data\n";
for(chrn=0; chrn<22; chrn++)
  {
  for( comp_it = compile_map[chrn].begin(); comp_it != compile_map[chrn].end(); ++comp_it )
    {  // Correcting and printing the corrected values
    for( unsigned snp_it=0; snp_it<12; snp_it++ ) if( (*comp_it).second.snp_pair == snp_pairs[snp_it] )
      for( unsigned sample=0; sample<num_rnaseq_samples; sample++ )
        if( AR_bias_supervec[ sample ][ snp_it ].median_AR > -1 )
          {  // Corrects the pvals and ARs using the SampleSnp median
          (*comp_it).second.sampleSnp_vec[ sample ].
            calc_pval ( AR_bias_supervec[ sample ][ snp_it ].median_AR );
          (*comp_it).second.sampleSnp_vec[ sample ].
            correct_AR( AR_bias_supervec[ sample ][ snp_it ].median_AR );
          }

    if( (*comp_it).second.usable_snp() )
      {
      (*comp_it).second.writeout( trim_out );
      for( unsigned sample=0; sample<num_rnaseq_samples; sample++ )
        if( (*comp_it).second.sampleSnp_vec[ sample ].good_het() )
          benj_vec[ sample ].push_back( (*comp_it).second.sampleSnp_vec[ sample ].pval );
      }
    }
  }

// Finally writeup the Benjamini-Hochberg thresholds for each sample
cout << "Calculating and writing the per sample Benjamini-Hochberg thresholds\n";
for( unsigned sample=0; sample < num_rnaseq_samples; sample++ )
  {
  vecsize = benj_vec[sample].size();
  sort( benj_vec[sample].begin(), benj_vec[sample].end() );
  for( unsigned i=1; i<vecsize; i++ )
    if( benj_vec[sample][i-1]*vecsize/(double)i > benj_fdr )
      { benj_out.file << " " << rnaseq_vec[sample] << "\t" << benj_vec[sample][i] << "\n"; break; }
  }


cout << "Done." << endl;

forbidden_in.tanca(); genoty_in.tanca();
trim_out.tanca(); fullinfo_out.tanca(); fullnoncl_out.tanca(); AR_bias_out.tanca();
}







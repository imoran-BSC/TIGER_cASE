/*
 * Allelic Imba: Select Signif
 *
 * This reads the Allelic_Imba_full/nonclonal_chrXX.txt files (twice), calculates the alignment
 * biases for all dinucleotide pairs in each sample, and corrects it in both AR and pval,
 * writing out only the nominally significantly imbalanced variants.
 *
 * This used to be the last part of filter_and_unify but was broken out to reduce mem footprint.
 *
 * Last modified: 13/11/2019
 */

#include "../includes.h"
#include "../class_fstream.h"
#include "../allelic_imba_common.h"  // Contains the functions to calculate imba, covg criteria

using namespace std;

bool debug = 0;  // DEBUG FLAG
const unsigned min_gtype_consist = 3;  // Min number of samples to check genotype consistency with RNAseq
const double benj_fdr = 0.01;  // 1% FDR Benjamini-Hochberg threshold per sample

unsigned num_samples;  // This substitutes num_rnaseq_samples and num_genoty_samples  !!

double binom_pval( const pair<int, int>& pr, const double p=0.5 );

/*
 *   CLASSES
 */

class class_count_stats
{  // Contains the ref/alt counts and vector of ARs
public:
  long int refc, altc;
  double raw_AR, mean_AR, median_AR;
  vector< double > AR_vec;

  class_count_stats() : refc(0), altc(0), raw_AR(-1), mean_AR(-1), median_AR(-1) {}
};

class class_sampleSnp
{  // Conte la info del SNP en aquesta mostra, la quantificacio d'allelic imbalance i el binom pval
public:
  string genoty;
  long int pose, refcount, altcount, errcount;
  double imba_val;

  class_sampleSnp() : genoty(""), pose(0), refcount(0), altcount(0), errcount(0), imba_val(0) {}

  bool bad_nclonal()
    { // Checks if there are enough clonal reads in this snp, returns false otherwise
    return( refcount + altcount < min_num_clon_reads ); }

  bool good_covg(double covg_mult=1)  // Checks if good coverage
    { return refcount + altcount >= (long int)(covg_mult * min_total_count); }

  bool good_het(double covg_mult=1)
    {  // If this sample is Het, with good coverage and AR
    return genoty.substr(0,3) == "Het" && good_covg(covg_mult) &&
      imba_val > imba_het_threshold && imba_val < 1-imba_het_threshold;
    }
  bool good_ref(double covg_mult=1)
    { return genoty == "Ref" && good_covg(covg_mult) && imba_val >= 1-imba_het_threshold; }
  bool good_alt(double covg_mult=1)
    { return genoty == "Alt" && good_covg(covg_mult) &&
        imba_val > -1 && imba_val <= imba_het_threshold; }

  void make_negative()
    { refcount = -1; altcount = -1; imba_val = -1; errcount = -1; }  // Destroy the clonal snp

  double calc_pval( double expected_ratio )
    {  // Calculates the binomial pvalue with the given expected allelic ratio
    double pval = -1;
    if( genoty.substr(0,3) == "Het" )  // Works for Het/Het1/Het2
      pval = binom_pval( pair<int,int>( refcount, altcount ), expected_ratio );
    return(pval);
    }

  void correct_AR( double expected_ratio )
    {  // Calculates the corrected ARs so that they are centered at 0.50 .
    double imba = -1;
    if( refcount+altcount > 0 )
      {
      imba = refcount/(double)(refcount+altcount);
      if( imba <= expected_ratio ) imba = imba*0.50/expected_ratio;
      else imba = 1-( (1-imba)*0.50/(1-expected_ratio) );
      }

    imba_val = imba;
    }
};

class class_compiledSnp
{  // Conte la info de totes les mostres per aquest SNP
public:
  string chr, ref, alt, rsname;
  long int pose;
  vector<class_sampleSnp> sampleSnp_vec;  // Vector of RNAseq'd samples SNP info

  class_compiledSnp() : chr(""), ref(""), alt(""), rsname(""), pose(0) {}

  void read_and_process( class_fstream& in )
    {  // Initialization
    class class_sampleSnp s;
    sampleSnp_vec.clear();
    sampleSnp_vec.resize( num_samples, s );  // Initialized to empty

    in.file >> chr >> pose >> rsname >> ref >> alt;
    for(unsigned k=0; k<num_samples; k++) {
      in.file >> s.refcount >> s.altcount >> s.errcount >> s.imba_val >> s.genoty;
      sampleSnp_vec[k] = s; }
    }

  unsigned get_snp_indx( vector< string >& snp_pairs )
    {  // Returns the di-nucleotide indx of this particular snp
    unsigned s;
    for( s=0; s<12; s++ )  // 12 dinucleotide permutations
      if( ref+alt == snp_pairs[s] )  // Separates SNPs by allelic pairs
        break;

    return s;
    }

  bool bad_sample_clonality( map<string, set<unsigned> >& badnc_rsids, unsigned& sample )
    {  // Bool to check if this particular rsid+sample combo is in the map
    map<string, set<unsigned> >::iterator map_it = badnc_rsids.find( rsname );
    if( map_it == badnc_rsids.end() ) return(false);  // SNP not in the list, so we are good
    else if( (*map_it).second.find( sample ) == (*map_it).second.end() ) return(false);  // SNP in the list, but this sample is good
    return(true);  // SNP in the list and this sample is bad
    }

  void negate_clonality( map<string, set<unsigned> >& badnc_rsids )
    {  // Puts -1s in the sampleSNPs prev detected as clonal
    map<string, set<unsigned> >::iterator map_it = badnc_rsids.find( rsname );

    for( unsigned sample=0; sample < num_samples; sample++)
      if( map_it != badnc_rsids.end() && (*map_it).second.find( sample ) != (*map_it).second.end() )
        sampleSnp_vec[ sample ].make_negative();  // Remove the clonal sampleSNPs
    }

  bool enough_usable_samples()
    {  // Checks if it has 3+ Het SNPs with min_num_reads+, reasonable ARs, and consistent genotypes
    bool good_genotype=0;
    unsigned num_ngt=0, num_good_hets=0;
    double refm=0, hetm=0, altm=0;
    vector<double> refvec, hetvec, altvec;  // Contain ALL Ref, Het and Alt samples AR

    for( unsigned k=0; k<num_samples; k++ )
      {
      if( sampleSnp_vec[k].genoty == "Ref" && sampleSnp_vec[k].imba_val > -1 )
        refvec.push_back( sampleSnp_vec[k].imba_val );
      else if( sampleSnp_vec[k].genoty.substr(0,3) == "Het" && sampleSnp_vec[k].imba_val > -1 )
        {
        if( sampleSnp_vec[k].good_het() ) num_good_hets++;
        hetvec.push_back( sampleSnp_vec[k].imba_val );
        }
      else if( sampleSnp_vec[k].genoty == "Alt" && sampleSnp_vec[k].imba_val > -1 )
        altvec.push_back( sampleSnp_vec[k].imba_val );
      else if( sampleSnp_vec[k].genoty == "" ) num_ngt++;
      }

    // Genotype consistency check by checking the median imbas versus reasonable ARs
    if( refvec.size() >= min_gtype_consist ) refm = funct_get_median( refvec ); else refm = -1;
    if( hetvec.size() >= min_gtype_consist ) hetm = funct_get_median( hetvec ); else hetm = -1;
    if( altvec.size() >= min_gtype_consist ) altm = funct_get_median( altvec ); else altm = -1;

    if( num_ngt != num_samples &&
        ( hetm > imba_het_threshold && hetm < 1-imba_het_threshold ) &&
        ( refm > 1-imba_het_threshold || refm < 0 ) &&  // To allow for non Ref samples
        ( altm < imba_het_threshold ) )
      good_genotype = 1;

    if( num_good_hets >= min_num_hets && good_genotype == 1 ) return 1;
    return 0;
    }

  void writeout( vector<double>& pval_vec, class_fstream& full_out )
    {  // Writes out (d'oh)
    string gty;

    full_out.file << chr << "\t" << pose << "\t" << rsname << "\t" << ref << "\t" << alt;
    for( unsigned k=0; k<num_samples; k++ )
      {  // Prints everything
      if( sampleSnp_vec[k].genoty == "" ) gty = "NGT";
      else gty = sampleSnp_vec[k].genoty;
      full_out.file << "\t" << sampleSnp_vec[k].refcount << "\t" << sampleSnp_vec[k].altcount << "\t"
        << sampleSnp_vec[k].errcount << "\t" << sampleSnp_vec[k].imba_val << "\t"
        << pval_vec[k] << "\t" << gty;
      }
    full_out.file << "\n";
    }
};

/*
 *   FUNCIONS SECUNDARIES
 */

void blob_of_text()
{  // A blob of text for when the wrong arguments are passed
cout << "\nThis program reads the full/nonclonal by-chr files produced by filter_and_unify, "
  << "calculates and corrects their dinucleotide biases, and writes the significant CASEs.\n\n"
  << "It requires 2 arguments:\n"
  << " 1) The path to where the 44 files (22 for each full/nonclonal) are located, with the "
  << "format Allelic_Imba_full/nonclonal_chrXX.txt, ending with a /"
  << "Ie /project/jferrer/imperial/Fastq_Data/input/\n"
  << " 2) The path to where the output files will be located, ending with a /"
  << "format Allelic_Imba_full/nonclonal_chrXX.txt"
  << "Ie /project/jferrer/imperial/Fastq_Data/output/\n"
  << "\nFor example, a valid call would be\n"
  << "$ ./select_signif /project/jferrer/imperial/Fastq_Data/input/ "
  << "/project/jferrer/imperial/Fastq_Data/output/\n"
  << "\nIf any of the above expectations are not met, the program will most probably "
  << "simply crash without much explanation or get stuck in an infinite loop.\n"
  << "\nThere are three output files:\n"
  << " - Allelic_Imba_trimmed.txt , which contains, for only those SNPs with usable Het values in 3+ samples, "
  << "the reference and alternate read counts, the corrected allelic ratio and p-value and their genotype.\n"
  << " - Allelic_Imba_AR_biases.txt , which contains the mean and median allelic ratio biases "
  << "per sample per dinucleotide.\n\n"
  << " - Allelic_Imba_Benjamini.txt , which contains the Benjamini-Hochberg significance threshold "
  << "per given input sample.\n\n"
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
int m=pr.first, n=pr.second, i=0, j=0;
double epsilon = 1e-7, b=0, limit = binom_prob( m, n, p ), cumul=0;

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

void funct_make_badnc_set( vector<class_fstream*>& invec, map<string, set<unsigned> >& badnc_rsids )
{  // Goes through the nonclonal files, and makes a set of the discarded snps by clonality
class_compiledSnp snp;

// Adds up all usable het ref and alt reads per sample and/or snp
for(unsigned chrn=0; chrn<22; chrn++)
  {
  snp.read_and_process( (*invec[chrn]) );
  while( !(*invec[chrn]).file.eof() )
    {
    for( unsigned sample=0; sample<num_samples; sample++ )
      if( snp.sampleSnp_vec[ sample ].bad_nclonal() )
        badnc_rsids[ snp.rsname ].insert(sample);  // We grab the bad sample for this SNP

    snp.read_and_process( (*invec[chrn]) );
    }
  }

cout << "Identified " << badnc_rsids.size() << " SNPs with 1+ clonal samples\n";
}

void funct_calc_AR_medians_and_writeout(
  vector<class_fstream*>& invec, vector< string >& snames_vec,
  vector< string >& snp_pairs, map<string, set<unsigned> >& badnc_rsids,
  vector< vector< class_count_stats > >& stats_vec,
  class_fstream& out )
{  // Calculates the per sample per dint allelic ratio biases and outputs them
unsigned snp_indx=0, vecsize=0, clonclounter=0;
long int refc=0, altc=0;
double m=0;
map<string, set<unsigned> >::iterator map_it;

class_compiledSnp snp;

set<string>::iterator set_it;
// Adds up all usable het ref and alt reads per sample and/or snp
// Per Sample per SNP calculate rawm mean and median ARs and printout

for(unsigned chrn=0; chrn<22; chrn++)
  {
  snp.read_and_process( (*invec[chrn]) );
  while( !(*invec[chrn]).file.eof() )
    {
    snp_indx = snp.get_snp_indx( snp_pairs );  // Select the correct dinucleotide index
    for( unsigned sample=0; sample<num_samples; sample++ )
      {
      if( snp.bad_sample_clonality( badnc_rsids, sample ) ) clonclounter++;
      else if( snp.sampleSnp_vec[ sample ].good_het() )  // and Het matches genoty+good_covg
        {
        refc = snp.sampleSnp_vec[ sample ].refcount;
        altc = snp.sampleSnp_vec[ sample ].altcount;

        // Compiles the number of ref and alt reads per sample per snp
        stats_vec[ sample ][ snp_indx ].refc += refc;
        stats_vec[ sample ][ snp_indx ].altc += altc;
        stats_vec[ sample ][ snp_indx ].
          AR_vec.push_back( refc/(double)(refc+altc) );
        }
      }

    snp.read_and_process( (*invec[chrn]) );
    }
  cout << "Done with chr" << chrn+1 << endl;
  }
cout << "Discarded a total of " << clonclounter << " snps due to clonality\n";

// Now for every sample+dnt combination, calculate and write the raw, mean and median AR
for( unsigned sample=0; sample<num_samples; sample++ )
  for( unsigned snp_it=0; snp_it<12; snp_it++ )
    {
    refc = stats_vec[sample][snp_it].refc;
    altc = stats_vec[sample][snp_it].altc;

    stats_vec[sample][snp_it].raw_AR = -1;
    if( refc+altc > 0 ) stats_vec[sample][snp_it].raw_AR = refc/(double)( refc+altc );

    m = 0; vecsize = stats_vec[sample][snp_it].AR_vec.size();
    stats_vec[sample][snp_it].mean_AR = -1;
    if(vecsize > 0) {
      for( unsigned k=0; k<vecsize; k++ ) m += stats_vec[sample][snp_it].AR_vec[k];
      stats_vec[sample][snp_it].mean_AR = m/(double)vecsize; }

    stats_vec[sample][snp_it].median_AR = funct_get_median( stats_vec[sample][snp_it].AR_vec );

    out.file << snames_vec[sample] << "_" << snp_pairs[snp_it] << "\t"
      << refc << "\t" << altc << "\t"  // Added as a sanity check?
      << stats_vec[sample][snp_it].raw_AR << "\t"
      << stats_vec[sample][snp_it].mean_AR << "\t" << stats_vec[sample][snp_it].median_AR << "\n";
    }
}

void funct_calc_correct_AR_pvals_and_writeout(
  vector<class_fstream*>& invec, vector< string >& snames_vec,
  vector< string >& snp_pairs, map<string, set<unsigned> >& badnc_rsids,
  vector< vector< class_count_stats > >& stats_vec,
  class_fstream& trim_out, class_fstream& benj_out )
{  // Reads the files again, corrects the AR biases, calculates the pvals and writeout
   // Now added the genoty/RNAseq consistency check
bool warn = false;
char header[65536];
unsigned snp_indx=0, vecsize=0;
double median_AR=0;

vector< unsigned > genoty_rnaseq_match(num_samples,0), genoty_rnaseq_total(num_samples,0);
vector< double > pval_vec;
vector< vector< double > > benj_vec(num_samples);
class_compiledSnp snp;


for(unsigned chrn=0; chrn<22; chrn++)
  {
  (*invec[chrn]).file.clear(); // clear bad state after eof
  (*invec[chrn]).file.seekg(0);
  (*invec[chrn]).file.getline(header, 65536);  // Process the header

  snp.read_and_process( (*invec[chrn]) );
  while( !(*invec[chrn]).file.eof() )
    {  // Correcting and printing the corrected values
    pval_vec.resize(num_samples, -1);

    // First we remove the samplesnps that were detected as clonal
    snp.negate_clonality( badnc_rsids );

    snp_indx = snp.get_snp_indx( snp_pairs );
    for( unsigned sample=0; sample<num_samples; sample++ )
      {  // Then correct the AR/pval values given the biases
      median_AR = stats_vec[ sample ][ snp_indx ].median_AR;  // Clonality was already removed
      if( median_AR > -1 )
        {  // Corrects the pvals and ARs using the SampleSnp median
        snp.sampleSnp_vec[ sample ].correct_AR( median_AR );
        pval_vec[sample] = snp.sampleSnp_vec[ sample ].calc_pval( median_AR );
        }
      }

    if( snp.enough_usable_samples() )  // Will print -1s in clonal samples
      {  // Finally if after correction the SNP is good, we writeout+benjamini
      snp.writeout( pval_vec, trim_out );
      for( unsigned sample=0; sample<num_samples; sample++ )
        {
        if( snp.sampleSnp_vec[ sample ].good_het() )
          benj_vec[ sample ].push_back( pval_vec[sample] );

        // We can now make a imba/genoty check per sample
        if( snp.sampleSnp_vec[sample].good_ref() ||
            snp.sampleSnp_vec[sample].good_het() ||
            snp.sampleSnp_vec[sample].good_alt() )
          { genoty_rnaseq_match[sample]++; genoty_rnaseq_total[sample]++; }
        else if( ( snp.sampleSnp_vec[sample].genoty == "Ref" ||
                   snp.sampleSnp_vec[sample].genoty.substr(0,3) == "Het" ||
                   snp.sampleSnp_vec[sample].genoty == "Alt" ) &&
                 snp.sampleSnp_vec[sample].good_covg() )  // Clonal samples have -1 coverage
          genoty_rnaseq_total[sample]++;  // If genotyped and good covg but no match
        }
      }

    snp.read_and_process( (*invec[chrn]) );
    }
  }

cout << "Calculating and writing the per sample Benjamini-Hochberg thresholds\n";
for( unsigned sample=0; sample < num_samples; sample++ )
  {
  vecsize = benj_vec[sample].size();
  sort( benj_vec[sample].begin(), benj_vec[sample].end() );
  for( unsigned i=1; i<vecsize; i++ )
    if( benj_vec[sample][i-1]*vecsize/(double)i > benj_fdr )
      { benj_out.file << " " << snames_vec[sample] << "\t" << benj_vec[sample][i] << "\n"; break; }
  }

// Sanity check comparing RNAseq and genotypes
for( unsigned sample=0; sample<num_samples; sample++ )
  if( genoty_rnaseq_match[sample]/(double)genoty_rnaseq_total[sample] < 0.90 )
    { warn = true; break; }
if( warn )
  {
  cout << "\nWARNING: Big discrepancy between RNA and genotype detected - "
    << "please make sure that the genotype file follows the same order than the "
    << "comma separated list of samples. The expected concordance values are >90%, "
    << "values around 50-60% indicate an unpaired sample. The values were:\n";
  for( unsigned sample=0; sample<num_samples; sample++ )
    cout << snames_vec[sample] << "\t"
      << genoty_rnaseq_match[sample] << "\t" << genoty_rnaseq_total[sample] << "\t"
      << genoty_rnaseq_match[sample]*100/(double)genoty_rnaseq_total[sample] << "%\n";
  cout << endl;
  }
}


/*   =====================
 *   === MAIN FUNCTION ===
 *   =====================
 */

int main( int argc, char* argv[] )
{
unsigned chrn=0;
char header[65536];
string inpath, outpath,
  cromosomes[22] =
  {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
   "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
   "chr20","chr21","chr22"},
  nucleotides[4]={"A","C","G","T"};

class_fstream AR_bias_out, benj_out, trim_out;
class_fstream* hi_inp;
vector<class_fstream*> full_invec, nc_invec;

vector< string > snp_pairs, snames_vec;
vector< vector< double > > benj_vec;
vector< vector< class_count_stats > > stats_vec;
map<string, set<unsigned> > badnc_rsids;

for(unsigned i=0; i<4; i++) for(unsigned j=0; j<4; j++) if(i != j)
  snp_pairs.push_back( nucleotides[i]+nucleotides[j] );  // Creates the 12 dinucleotide pairs

if(!debug) cout << "Starting SELECT SIGNIF" << endl;
else cout << "Starting SELECT SIGNIF -- in DEBUG MODE --" << endl;
{  // Obre totes les entrades i sortides
if( debug )
  {
  // This is both the input an
  inpath = "/home/bscuser/project_data/t2dsystems/case/pipeline_tests/test_hiega/";
//  inpath = "/home/ender/projects_data/test/";
  outpath = inpath;
  }
else if( !debug && argc != 3 ) { blob_of_text(); cout << "Terminating.\n"; return 1; }
else if( !debug && argc == 3 )
  {  // Arguments from command line
  inpath  = argv[1];
  outpath = argv[2];

  cout << "Given paramters were:\n";
  for(unsigned u=1; u<unsigned(argc); u++) cout << u << ") " << argv[u] << "\n";
  cout << "\n";
  }

// Load the full/nonclonal by-chr files in two separated pointer vectors
for( chrn=0; chrn<22; chrn++ )
  {  // Pointer vecs to all quantification files, standard and non-clonal
  try { hi_inp = new class_fstream; }  // Using NEW and DELETE so that the file does not go out of scope
  catch (bad_alloc& ba) { cout << "bad_alloc caught: " << ba.what() << "\n"; }

  if( !hi_inp->obre( inpath + "Allelic_Imba_full_" + cromosomes[chrn] + ".txt",
      "in", "full_" + cromosomes[chrn] ) )
    full_invec.push_back( hi_inp );
  else return 1;

  // Now non-clonal
  try { hi_inp = new class_fstream; }  // Using NEW and DELETE so that the file does not go out of scope
  catch (bad_alloc& ba) { cout << "bad_alloc caught: " << ba.what() << "\n"; }

  if( !hi_inp->obre( inpath + "Allelic_Imba_nonclonal_" + cromosomes[chrn] + ".txt",
      "in", "nonclonal_" + cromosomes[chrn] ) )
    nc_invec.push_back( hi_inp );
  else return 1;
  }

// We need to grab the name of all samples from the header of the 1st file
// This ASSUMES that each file contains the info from the same samples, in the same order
// (which is what should happen in a clear run of the pipeline)
string strheader, tmp;
int field=0; size_t p1=0, p2=0, p3=0;

hi_inp = full_invec[0];  // Grab the 1st file from the pointer vec
(*hi_inp).file.getline(header, 65536);
strheader = string(header);

while( ( p2 = strheader.find( "\t", p1 ) ) != string::npos )
  { if( field > 4 && (field-5) % 5 == 0 ) {
      tmp = strheader.substr( p1, p2-p1 ); p3 = tmp.find("_refc");
      snames_vec.push_back( tmp.substr(0,p3) ); }
    field++; p1=p2+1; }

num_samples = snames_vec.size();
cout << "We have read " << num_samples << " samples with data to process\n";

// And a few vector sizes that depend on it
stats_vec.resize(
  num_samples, vector< class_count_stats >(12) );  // 12 dinucleotide combinations
benj_vec.resize( num_samples );

// Now open the output files, write headers
if( !AR_bias_out.obre( outpath + "Allelic_Imba_AR_biases.txt", "out" ) &&  // Per-sample per-dinucleotide allelic ratio biases, for posterior correction
    !benj_out.obre   ( outpath + "Allelic_Imba_Benjamini.txt", "out" ) &&  // Per-sample Benjamini Hockberg thresholds
    !trim_out.obre   ( outpath + "Allelic_Imba_trimmed.txt", "out" ) )  // Trimmed results to only 3+ valid heterozgous SNPs
  {
  stringstream ss;
  // Read the 1st line of all inputs, headers
  for(unsigned k=0; k<22; k++)
    {
    if(k > 0) {  // As we have already read the 1st file's header
      hi_inp = full_invec[k];  // Grab the 1st file from the pointer vec
      (*hi_inp).file.getline(header, 65536); }

    hi_inp = nc_invec[k];  // Also process all nonclonal files
    (*hi_inp).file.getline(header, 65536);
    }

  // Write output headers
  AR_bias_out.file << "#Name\tTotalRef\tTotalAlt\tTotalAR\tMeanAR\tMedianAR\n";

  ss << "#chr\tposition\trsname\tref\talt";
  for( unsigned k=0; k<num_samples; k++ )
    ss << "\t" << snames_vec[k] << "_refc\t" << snames_vec[k] << "_altc\t"
      << snames_vec[k] << "_errc\t" << snames_vec[k] << "_imba\t"
      << snames_vec[k] << "_pval\t" << snames_vec[k] << "_genoty";
  ss << "\n";
  trim_out.file << ss.str();  // We now add the _pval field to this file
  }
else return 1;
}

// First go through the nonclonal files, and make a set containing the snps that pass the clonality filters
cout << "Discarding SNPs by clonality" << endl;
funct_make_badnc_set( nc_invec, badnc_rsids );


// Now let's read all full files, and calculate the total, mean and median ARs
// for each sample's dinucleotide pairs
cout << "Calculating and writing the raw AR biases" << endl;
funct_calc_AR_medians_and_writeout(
  full_invec, snames_vec, snp_pairs, badnc_rsids, stats_vec, AR_bias_out );


// Now we need to apply the corrections to the observed data, calculate the (corrected) pvals, and writeout
cout << "Correcting the bias and writing the trimmed SNP data" << endl;
funct_calc_correct_AR_pvals_and_writeout(
  full_invec, snames_vec, snp_pairs, badnc_rsids, stats_vec, trim_out, benj_out );



cout << "All done!" << endl;

AR_bias_out.tanca(); benj_out.tanca(); trim_out.tanca();

return(0);
}

























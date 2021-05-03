/*
 *  REGULOME COMMON -BODY-
 *
 *  Contains the current SNP coverage criteria (ie >10 reads min allele) and a function to test it,
 *  as well as current definition of Allelic Imbalance Likelihood and a function to calculate it,
 *  and a function to calculate and return the binomial pvalue distribution of the selected SNPs
 *
 */


#include "includes.h"
#include "class_fstream.h"
#include "allelic_imba_common.h"

using namespace std;

unsigned num_rnaseq_samples,
  num_genoty_samples;  // Exported from the header so it is visible everywhere else

/*
 *   DEFINITION OF CLASS FUNCTIONS
 */

void class_snp::read( class_fstream& in )
{  // Reads the SNP info from all samples
num_hets=0;
sample_vec.clear(); sample_vec.resize( num_rnaseq_samples );

in.file >> chr >> pose >> rsname >> ref >> alt;

// class_snp needs to be INITIALISED or else sample_vec is an empty vector!!
for( unsigned k=0; k<num_rnaseq_samples; k++ )
  {
  in.file >> sample_vec[k].refc >> sample_vec[k].altc >> sample_vec[k].errc
    >> sample_vec[k].imba >> sample_vec[k].pval >> sample_vec[k].genoty;
  if( sample_vec[k].genoty.substr(0,3) == "Het" ) num_hets++;
  if( sample_vec[k].pval == 0 ) sample_vec[k].pval = 1e-295;  // For the few zero cases
  if( sample_vec[k].refc+sample_vec[k].altc > 0 &&
      sample_vec[k].errc/(sample_vec[k].refc+sample_vec[k].altc) < err_read_pctg )
    sample_vec[k].totalcount = sample_vec[k].refc+sample_vec[k].altc;
  else sample_vec[k].totalcount = -1;
  }
hets_imba_pvec.clear(); hets_pval_pvec.clear(); hets_readc_pvec.clear();
}

void class_snp::calc_median_coverage()
{  // Calculates the median coverage of the viable reporter het SNPS
vector< double > count_vec;

for( unsigned k=0; k<num_rnaseq_samples; k++ )
  if( funct_coverage_test( *this , k ) )
    count_vec.push_back( sample_vec[k].refc + sample_vec[k].altc );

if( count_vec.size() >= min_num_hets )
  median_coverage = funct_get_median( count_vec );  // Median amount of read coverage of the SNP
else median_coverage = -1;
}

void class_snp::load_median_info( map<double, pair<int,double> >& median_map, vector<double>& median_list )
{  // Loads the median map and list
map<double, pair<int,double> >::iterator map_it;

calc_median_coverage();

if(num_hets >= min_num_hets)
  {
  map_it = median_map.find( median_coverage );
  if( map_it != median_map.end() )
    { ( (*map_it).second ).first++; ( (*map_it).second ).second++; } // Increments the count of the map
  else median_map[ median_coverage ] = pair<int,double>(1,1);  // Iniciates the median map count

  if( median_coverage > 0 )
    median_list.push_back( median_coverage );
  }
}

void class_snp::load_shuffle_vec( vector<class_sampledata>& shuffle_vec )
{  // Uploads the data to the shuffle vector
for( unsigned k=0; k<num_rnaseq_samples; k++ )
  if( sample_vec[k].genoty.substr(0,3) == "Het" )
    shuffle_vec.push_back( sample_vec[k] );
}

bool class_snp::operator< ( const class_snp& s ) const
{ return pose < s.pose; }

void class_genotype::read( class_fstream& in )
{  // Llegeix la info de la entrada i simplifica el genotipatge
string dummy;
vector<string> nt_vec( num_genoty_samples );
genoty_vec.clear(); genoty_vec.resize( num_genoty_samples );

in.file >> chr >> pose >> rsname >> ref >> alt >> dummy;

for( unsigned k=0; k<num_genoty_samples; k++ )
  {  // Converteix el genotipatge de NTs a Ref/Het/Alt/NAN
  in.file >> nt_vec[k];

  if( nt_vec[k][0] == ref[0] && nt_vec[k][1] == ref[0] ) genoty_vec[k] = "Ref";
  else if ( ( nt_vec[k][0] == ref[0] && nt_vec[k][1] == alt[0] ) ||
            ( nt_vec[k][0] == alt[0] && nt_vec[k][1] == ref[0] ) ) genoty_vec[k] = "Het";
  else if( nt_vec[k][0] == alt[0] && nt_vec[k][1] == alt[0] ) genoty_vec[k] = "Alt";
  else genoty_vec[k] = "NAN";
  }
}

void class_genotype::read_vcf( class_fstream& in )
{  // Reads and processes phased or unphased genotype data in VCF format
size_t p1;
string cchr, dummy, genoty;

genoty_vec.resize( num_genoty_samples );

in.file >> cchr >> pose >> rsname >> ref >> alt >> dummy >> dummy >> dummy >> dummy;
chr = "chr"+cchr;  // Since vcf format only contains numeric chromosomes

for( unsigned k=0; k<num_genoty_samples; k++ )
  {  // Read VCF genotype and convert it to Ref/Het/Alt/NAN
  in.file >> dummy;

  // Includes phased VCF format!
  p1 = dummy.find( ":" );
  genoty = dummy.substr( 0, p1 );

  genoty_vec[k] = "NAN";
  if     ( genoty == "0|0" || genoty == "0/0" ) genoty_vec[k] = "Ref";
  else if( genoty == "0|1" ) genoty_vec[k] = "Het1";
  else if( genoty == "1|0" ) genoty_vec[k] = "Het2";
  else if( genoty == "0/1" ) genoty_vec[k] = "Het";
  else if( genoty == "1|1" || genoty == "1/1" ) genoty_vec[k] = "Alt";
  }
}

bool class_genotype::operator< ( const class_genotype& g ) const
{ return pose < g.pose; }

/*
 *   FUNCIONS SECUNDARIES
 */

double funct_get_median( vector<double>& vec )
{  // Returns the median of a given vector
unsigned cnt=0, sz = vec.size();
double rt=0;
vector<double> v;
vector<double>::iterator it;

if( sz == 0 ) return -1;
else if( sz == 1 ) return vec[0];
else if( sz == 2 ) return (vec[0]+vec[1])/2;

//for( unsigned k=0; k<sz; k++ ) v.push_back( vec[k] );
v = vec;  // Copy since we will sort it and that would fuck up the original
sort( v.begin(), v.end() );

it = v.begin();
while( cnt < sz/(double)2-0.51 ) { cnt++; ++it; }
rt = (*it);  // Is the middle value of odd vectors, or middle+0.5 for even
if( sz%2 == 0 ) { --it; rt += (*it); rt /= 2; }  // Mean with middle-0.5 for even vectors

return rt;
}

bool funct_coverage_test( const class_snp& snp, const unsigned& indx )
{  // Contains the current "good enough heterozygous SNP" for bool testing
return snp.sample_vec[ indx ].genoty.substr(0,3) == "Het" &&    // If it is a RNA Het
  snp.sample_vec[ indx ].totalcount >= (int)min_total_count &&
  snp.sample_vec[ indx ].imba > imba_het_threshold &&
  snp.sample_vec[ indx ].imba < 1-imba_het_threshold ;    // To avoid absoulte 0 or 1s (artifacts?)
//  snp.sample_vec[ indx ].refc > min_allele_count &&
//  snp.sample_vec[ indx ].altc > min_allele_count ;    // To avoid absoulte 0 or 1s (artifacts?)
}

double funct_calc_fisher( vector<double*>& pv )
{  // Calculates the Fisher's X^2 statistics of the given pvals
unsigned vecsize = pv.size();
double xsq=0;

for( unsigned u=0; u<vecsize; u++ )
  {
  if( *(pv[u]) <= 0 ) return -1;  // Which doesnt happen for pvals since at worst they are 1
  if( *(pv[u]) < min_pval )
    xsq += log( min_pval ); // Necessary to avoid 1 result shifting it all
  else xsq += log( *(pv[u]) );  // Natural Log as per the formula
  }

return -2*xsq;
}

double funct_qnorm( const double& pv )
{  // Given a pvalue, returns the qnorm value from the table, after boundary checks
unsigned k;
double z;

if( pv < 1e-15 ) return qnorm_vals[qnorm_length-1];

for( k=1; k<qnorm_length; k++ )
  if( k/100.0 > -log10(pv) ) { z = qnorm_vals[k-1]; break; }  // pow( 10.0, -(k/100.0) ) <= pv
if( k == qnorm_length ) z = qnorm_vals[qnorm_length-1];  // Last value catchall, reduntant safety

return z;
}

double funct_calc_zvalue( vector<double*>& pv, vector<double*>& ratios,
  vector<int*>& read_counts, int with_direction )  // with_direction defaults to 1 in declaration
{  // Calculates the Stouffer's Z-score, with or without direction
unsigned good=0, vecsize = pv.size();
double d, up=0, down=0;

vector<double> z(vecsize);

for( unsigned i=0; i<vecsize; i++ )
  {
  if( *(pv[i]) <= 0 || *(ratios[i]) <=0 || *(read_counts[i]) <= 0 ) continue;  // Skip this one - should not happen

  z[i] = funct_qnorm( *(pv[i])/2.0 );  // Fetches the right Z-scores from the R table, z=qnorm(1-p/2)

  d=1;
//  if( with_direction && ( *(ratios[i]) < 0.5 || ( *(ratios[i]) == 0.5 && rand()/(double)RAND_MAX < 0.5 ) ) )
//    d = -1;  // Randomly assings a sign if AR exactly 0.5
  if( with_direction && *(ratios[i]) < 0.5 ) d = -1;
  if( *(ratios[i]) == 0.5 ) d = 0;  // Remove exactly 0.5 AR
  up += z[i]*sqrt( *(read_counts[i]) )*d;  // as per paper Z = sum(z_i * w_i * d_i )/sqrt( sum(w_i^2) )
  down += *(read_counts[i]);
  good++;
  }

if( down > 0 && good >= min_num_hets ) return up/sqrt( down );
else return -100;
}







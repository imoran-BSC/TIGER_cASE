/*
 *  PERMUTED REPORTER IMBALANCE
 *
 *  This script reads an Allelic_Imba_trimmed.txt file,
 *  calculates the empiric Zscore distribution, then randomises
 *  the Het read counts 1000 times within 5 bins (one for SNPs with
 *  a median coverage of zero, and four more containing 20% of the
 *  reminding data) and creates an expected null Zscore distribution.
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

bool debug = 0;  // DEBUG FLAG
const unsigned num_median_bins = 1+4,  // 2 for >0 and <0 bins, >2 for equally spaced positive bins
  num_quartiles = 5000;  // Number of points for the qqplot
const double fdr_levels[3] = { 2.5, 5, 10 };  // FDR values for which we will calculate and report the AIL value

unsigned num_permutations;  // Init in the config file

vector<double> median_bins;


/*
 *   CLASSES
 */

// class class_sampledata  // Defined in regulome_common.h: Contains refc, altc, errc, imba, pval, genoty
// class class_snp         // Defined in regulome_common.h: Contains all the SNP data + functions

/*
 *   FUNCIONS SECUNDARIES
 */

void funct_cumul_and_normalize( map<double, pair<int,double> >& data_map )
{  // Calculates the normalized cumulative map of a histogram
long int cumul=0;
map<double, int>::iterator hist_it;
map<double, pair<int,double> >::iterator data_it;

for( data_it=data_map.begin(); data_it!=data_map.end(); ++data_it )
  { ( (*data_it).second ).second += cumul;
    cumul = ( (*data_it).second ).second; }

for( data_it=data_map.begin(); data_it!=data_map.end(); ++data_it )
  ( (*data_it).second ).second = ( (*data_it).second ).second/cumul;  // Normalizing to 1
}

void funct_calc_median_bins_and_writeout( vector<double>& median_vec,
  map<double, pair<int,double> >& median_map, class_fstream& median_out )
{  // Calculate the bin median limits depending on the number of requested bins and writeout median hist
unsigned median_size, uit;
vector<double>::iterator median_it;
map<double, pair<int,double> >::iterator map_it;

if( num_median_bins > 2 )
  {
  sort( median_vec.begin(), median_vec.end() ); median_size = median_vec.size();
  median_it = median_vec.begin(); uit=0;
  median_bins.push_back( 0 );
  for( unsigned k=1; k<(num_median_bins-1); k++ )
    { for( ; median_it != median_vec.end() && uit < k*median_size/(double)( num_median_bins-1 ); uit++ )
        ++median_it;
      median_bins.push_back( *median_it ); }
  }
else if( num_median_bins == 2 ) median_bins.push_back( 0 );  // 2 bins are always <0 and >0
median_bins.push_back( 25000 );  // Last number is always way greater than possible

funct_cumul_and_normalize( median_map );  // Make the normalized cumulative

// Writeout the ordered median coverage map
median_out.file << "#Median breaks were";
for( unsigned u=0; u<median_bins.size(); u++ ) median_out.file << " " << median_bins[u];
median_out.file << "\n";
for( map_it = median_map.begin(); map_it != median_map.end(); ++map_it )
  median_out.file << (*map_it).first << "\t" << ( (*map_it).second ).first << "\t"
    << ( (*map_it).second ).second << "\n";
median_out.file.flush();
}

void funct_make_hist_map( map< double, pair<int,double> >& hist_map )
{  // Initiates the map with N bins within the given X-Y range
int lower_lim = -25, upper_lim =0, num_steps = 10000;  // Logarithmic10 limits to the map

for( int j=0; j<=num_steps; j++ )  // Iterate from -30 to 0 in 10000 parts
  hist_map[ pow( 10, lower_lim + j*( upper_lim-lower_lim )/(double)num_steps ) ] = pair<int,double>(0,0);
}

void funct_put_back_values( vector< list<class_snp> >& snp_vec,
  vector< vector<class_sampledata> >& shuffle_vecs,
  vector<unsigned>& shuffle_sizes )
{  // Fills the genotyped Hets of the snp list with the randomized values from the vector
int chrn=0;
list<class_snp>::iterator snp_it;
vector<unsigned> shuffle_it( num_median_bins, 0 );  // Initiates an interator for each bin at zero

for( chrn=0; chrn<22; chrn++ )
  for( snp_it = snp_vec[chrn].begin(); snp_it != snp_vec[chrn].end(); ++snp_it )
    if( (*snp_it).num_hets >= min_num_hets )
      for( unsigned b=0; b<num_median_bins; b++ )
        if( (*snp_it).median_coverage <= median_bins[b] )  // Selects the right shuffle vec
          {
          for( unsigned k=0; k<num_rnaseq_samples; k++ )
            if( shuffle_it[b] < shuffle_sizes[b] && (*snp_it).sample_vec[k].genoty.substr(0,3) == "Het" )
              {  // Enters new values for the genotyped Hets of this coverage bin
              (*snp_it).sample_vec[k] = shuffle_vecs[b][ shuffle_it[b] ];
              shuffle_it[b]++;
              }
          break;
          }
}

void funct_calc_qqplot_and_writeout( vector<double>& rnd_list, vector<double>& obs_list,
  class_fstream& rnd_distr_out, class_fstream& qqplot_out, int type )
{  // Writeout the qqplot according to the given number of quantiles (points in the graph)
unsigned choose, num_bins, max_val=50;  // MAX expected QQplot value
long unsigned rnd_norm=0, obs_norm=0, rnd_counter, obs_counter, rnd_size, obs_size;
double interval = 0.25;
vector<double> tmp_list, rnd_negative_list, obs_negative_list;
vector<double>::iterator tmp_it, rnd_it, obs_it;
map<double, pair<double, double> > distr_map;
map<double, pair<double, double> >::iterator map_it;

sort( rnd_list.begin(), rnd_list.end() ); sort( obs_list.begin(), obs_list.end() );

num_bins = max_val/interval;

// Makes a map of the distribution of the randomized and obs data
if( type == 1 )
  for(unsigned b=0; b<num_bins+1; b++)
    distr_map[ interval/(double)2+b*interval ] = pair<double,double>(0,0);
else if( type == 2 )  // Also create the negative part of the map
  for(unsigned b=0; b<num_bins*2+1; b++)
    distr_map[ (-1)*(double)max_val+interval/(double)2+b*interval ] = pair<double,double>(0,0);

// We can now randomly choose upper and lower bound syntaxes to get rid of boundary bias
for( rnd_it=rnd_list.begin(); rnd_it!=rnd_list.end(); ++rnd_it )
  {  // Populates the random map
  choose = rand() % 2;  // Either zero or one
  if( choose == 0 )
    {  // Use lower bound syntax
    map_it = distr_map.lower_bound( *rnd_it );
    if( map_it != distr_map.end() ) { (*map_it).second.first++; rnd_norm++; }
    }
  else if( choose == 1 )
    {  // Use upper bound syntax
    map_it = distr_map.upper_bound( *rnd_it );
    if( map_it != distr_map.end() ) { (*map_it).second.first++; rnd_norm++; }
    }
  }
for( obs_it=obs_list.begin(); obs_it!=obs_list.end(); ++obs_it )
  {  // Populates the obs map
  choose = rand() % 2;  // Either zero or one
  if( choose == 0 )
    {  // Use lower bound syntax
    map_it = distr_map.lower_bound( *obs_it );
    if( map_it != distr_map.end() ) { (*map_it).second.second++; obs_norm++; }
    }
  else if( choose == 1 )
    {  // Use upper bound syntax
    map_it = distr_map.upper_bound( *obs_it );
    if( map_it != distr_map.end() ) { (*map_it).second.second++; obs_norm++; }
    }
  }

// Normalization and printout
for( map_it = distr_map.begin(); map_it != distr_map.end(); ++map_it )
  rnd_distr_out.file << (*map_it).first-interval/(double)2 << "\t"
    << (*map_it).second.second/(double)(obs_norm*interval) << "\t"
    << (*map_it).second.first/(double)(rnd_norm*interval) << "\n";

// Sorts the rnd lists and writes out the quantiles
if( type == 2 )
  {  // Separate the negative results in different lists and make them positive, then rerun qqplot
  tmp_list.clear();
  for( obs_it = obs_list.begin(); obs_it != obs_list.end(); ++obs_it )
    {
    if( *obs_it >0 ) tmp_list.push_back( *obs_it );
    else if( *obs_it <0 ) obs_negative_list.push_back( (-1)*(*obs_it) );
    }  // Leaving the zeroes alone to remove bias
  sort( obs_negative_list.begin(), obs_negative_list.end() );
  obs_list.clear(); obs_list = tmp_list;

  tmp_list.clear();
  for( rnd_it = rnd_list.begin(); rnd_it != rnd_list.end(); ++rnd_it )
    {
    if( *rnd_it >0 ) tmp_list.push_back( *rnd_it );
    else if( *rnd_it <0 ) rnd_negative_list.push_back( (-1)*(*rnd_it) );
    }
  sort( rnd_negative_list.begin(), rnd_negative_list.end() );
  rnd_list.clear(); rnd_list = tmp_list;
  }
obs_it = obs_list.begin(); rnd_it = rnd_list.begin();
rnd_counter=0; obs_counter=0; rnd_size = rnd_list.size(); obs_size = obs_list.size();
for( unsigned q=1; q<=num_quartiles; q++ )
  {  // Iterates the sorted lists, and prints the rounded quartiles
  while( rnd_it != rnd_list.end() &&
         rnd_counter < (long unsigned)(rnd_size*q/(double)num_quartiles ) )
    { rnd_counter++; ++rnd_it; }
  while( obs_it != obs_list.end() &&
         obs_counter < (long unsigned)(obs_size*q/(double)num_quartiles ) )
    { obs_counter++; ++obs_it; }

  if( rnd_it == rnd_list.end() ) --rnd_it;
  if( obs_it == obs_list.end() ) --obs_it;

  if( *rnd_it != 1 && *obs_it != 1 )
    qqplot_out.file << *rnd_it << "\t" << *obs_it << "\n";
  // To fix the returning negative zeroes (is not 1 exactly but diff underflows)
  else if( *rnd_it == 1 && *obs_it != 1 )
    qqplot_out.file << "0\t" << *obs_it << "\n";
  else if( *rnd_it != 1 && *obs_it == 1)
    qqplot_out.file << *rnd_it << "\t0\n";
  else if( *rnd_it == 1 && *obs_it == 1)
    qqplot_out.file << "0\t0\n";
  }
if( type == 2 )
  {  // We redo the calculation for the negative vectors, with inverted looping
  obs_it = obs_negative_list.begin(); rnd_it = rnd_negative_list.begin();
  rnd_counter=0; obs_counter=0; rnd_size = rnd_negative_list.size(); obs_size = obs_negative_list.size();
  for( unsigned q=1; q<=num_quartiles; q++ )
    {  // Iterates the sorted lists, and prints the rounded quartiles
    while( rnd_it != rnd_negative_list.end() &&
           rnd_counter < (long unsigned)(rnd_size*q/(double)num_quartiles ) )
      { rnd_counter++; ++rnd_it; }
    while( obs_it != obs_negative_list.end() &&
           obs_counter < (long unsigned)(obs_size*q/(double)num_quartiles ) )
      { obs_counter++; ++obs_it; }

    if( rnd_it == rnd_negative_list.end() ) --rnd_it;
    if( obs_it == obs_negative_list.end() ) --obs_it;

    if( *rnd_it != 1 && *obs_it != 1 )
      qqplot_out.file << *rnd_it << "\t" << *obs_it << "\n";
    // To fix the returning negative zeroes (is not 1 exactly but diff underflows)
    else if( *rnd_it == 1 && *obs_it != 1 )
      qqplot_out.file << "0\t" << *obs_it << "\n";
    else if( *rnd_it != 1 && *obs_it == 1)
      qqplot_out.file << *rnd_it << "\t0\n";
    else if( *rnd_it == 1 && *obs_it == 1)
      qqplot_out.file << "0\t0\n";
    }
  }
}

/*   =====================
 *   === MAIN FUNCTION ===
 *   =====================
 */

int main( int argc, char* argv[] )
{
unsigned chrn=0, num_covg_hets=0;
char header[65536];
double zvalue;
string trimmed_file, out_path, strheader, chrt="", sys, cromosomes[22] =
  {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
   "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
   "chr20","chr21","chr22"};  // ChrX and Y are out

class_snp snp;
map<double, pair<int,double> > median_map;
vector<unsigned> shuffle_sizes( num_median_bins );
vector<double> median_vec, obs_zval_vec, rnd_zval_vec;
vector<int*> hets_readc_pvec;  // To store the randomized data pointers
vector<double*> hets_imba_pvec, hets_pval_pvec;
vector<string> hi_vec;
vector< list<class_snp> > snp_vec(22);
list<class_snp>::iterator snp_it;
vector< vector<class_sampledata> > shuffle_vecs(num_median_bins);  // One shuffle vec per median bin

class_fstream config_in, data_in, median_out, rnd_zval_out, qqplot_zval_out;

srand( unsigned( time( NULL ) ) );  // Random seed initialization

if(!debug) cout << endl << "Starting ALLELIC IMBA RANDOMIZER\n" << endl;
else cout << endl << "Starting ALLELIC IMBA RANDOMIZER -- in DEBUG MODE --\n" << endl;
{  // Open all inputs and outputs
unsigned num_arguments = 3;
if(!debug)
  {  // Read input arguments from config file
  if( !config_in.open( "config.ini", "in", "config" ) )
    {
    string config_argv[num_arguments+1], dummy;
    cout << "Read from config file:\n";
    for(unsigned k=1; k<num_arguments+1; k++){
      config_in.file >> dummy >> config_argv[k];
      cout << dummy <<"\t"<< config_argv[k] << "\n";}

    trimmed_file = config_argv[1];
    num_permutations = atoi( config_argv[2].c_str() );
    out_path = config_argv[3];
    }
  else { cout << "Error reading the config file, terminating!"; return 1; }
  config_in.close();
  }
else if(debug && argc == (int)num_arguments+1)
  {  // In debug mode, use command line arguments instead
    trimmed_file = argv[1];
    num_permutations = atoi( argv[2] );
    out_path = argv[3];
  }
else if(debug && argc != (int)num_arguments+1)
  {cout << "Wrong number of arguments. Terminating!"; while( !cin.get() ); return 1;}

if( !data_in.open( trimmed_file, "in", "data" ) &&

    !rnd_zval_out.open   ( out_path + "Allelic_Imba_rnd_Zscore.txt", "out" ) &&
    !qqplot_zval_out.open( out_path + "Allelic_Imba_QQplot_Zscore.txt", "out" ) &&
    !median_out.open     ( out_path + "Allelic_Imba_Cumul_median_covg.txt", "out" ) )
  {
  // Get the number and name of samples from the header of Allelic_Imba_trimmed.txt
  size_t p1=0, p2=0, p3=0;
  unsigned field=0;
  data_in.file.getline(header, 65536); strheader = string(header);
  while( ( p2 = strheader.find( "\t", p1 ) ) != string::npos )
    { field++; if( field % 6 == 0 ) { p3 = strheader.find( "_refc", p1 );
        hi_vec.push_back( strheader.substr( p1, p3-p1 ) ); } p1=p2+1; }

  num_rnaseq_samples = hi_vec.size();  // GLOBAL DEFINITION OF num_rnaseq_samples
  cout << "We have " << num_rnaseq_samples << " samples\n";
  }
else { return 1; }
}

cout << "Reading SNP data\n";
chrt=""; snp.read( data_in );
while( !data_in.file.eof() )
  {  // Reads all the SNP data
  if( chrt != snp.chr )
    {
    chrt=snp.chr;
    for(chrn=0; chrn<22; chrn++)
      if( cromosomes[chrn] == chrt ) break;
    }

  if( chrn<22 )
    {
    snp.load_median_info( median_map, median_vec );
    if( snp.num_hets >= min_num_hets ) snp_vec[chrn].push_back(snp);
    }

  snp.read( data_in );
  }

for(chrn=0;chrn<22;chrn++)
  for( snp_it=snp_vec[chrn].begin(); snp_it!=snp_vec[chrn].end(); ++snp_it )
    {  // Calculates the Zval and pushes the > -1 to the obs list for posterior shuffling
    for( unsigned k=0; k<num_rnaseq_samples; k++ )
      if( funct_coverage_test( (*snp_it), k ) )  // If it abides current read depth standard
        { (*snp_it).hets_imba_pvec.push_back( &(*snp_it).sample_vec[k].imba );
          (*snp_it).hets_pval_pvec.push_back( &(*snp_it).sample_vec[k].pval );
          (*snp_it).hets_readc_pvec.push_back( &(*snp_it).sample_vec[k].totalcount ); }

    (*snp_it).num_covg_hets = (*snp_it).hets_imba_pvec.size();

    if( (*snp_it).num_covg_hets >= min_num_hets )
      {  // Calculate Zval and if consistent or not
      (*snp_it).zvalue = funct_calc_zvalue( (*snp_it).hets_pval_pvec,
        (*snp_it).hets_imba_pvec, (*snp_it).hets_readc_pvec );
      }

    if( (*snp_it).zvalue > -100 ) obs_zval_vec.push_back( (*snp_it).zvalue );
    }

// Calculate the bin median limits depending on the number of requested bins and writeout median hist
funct_calc_median_bins_and_writeout( median_vec, median_map, median_out );

// Now separate the SNPs in the different shuffle vecs
for(chrn=0; chrn<22; chrn++)
  for( snp_it=snp_vec[chrn].begin(); snp_it!=snp_vec[chrn].end(); ++snp_it )
    for( unsigned b=0; b<num_median_bins; b++ )
      if( (*snp_it).median_coverage <= median_bins[b] )
        { (*snp_it).load_shuffle_vec( shuffle_vecs[b] ); break; }  // Loads the data to the right shuffle vec

cout << "\nMedian breaks at";
for( unsigned b=0; b<num_median_bins; b++ ) cout << " " << median_bins[b];
cout << " coverage\nPerforming " << num_permutations << " permutations.\n";


for( unsigned b=0; b<num_median_bins; b++ )
  shuffle_sizes[b] = shuffle_vecs[b].size();

for( unsigned k=0; k<num_permutations; k++)
  {  // Performing the permutations and imbalance calculations
  if( (k+1)%10==0 ) cout << "Permutation " << k+1 << " of " << num_permutations << "...\n";
  for( unsigned b=0; b<num_median_bins; b++ )
    random_shuffle( shuffle_vecs[b].begin(), shuffle_vecs[b].end() );  // Shuffles the Het info vectors

  // Put back the shuffled values to their mean read coverage corresponding genotyped Hets
  funct_put_back_values( snp_vec, shuffle_vecs, shuffle_sizes );

  // And calculate their Zval
  for(chrn=0; chrn<22; chrn++)
    for( snp_it=snp_vec[chrn].begin(); snp_it!=snp_vec[chrn].end(); ++snp_it )
      {
      hets_imba_pvec.clear(); hets_pval_pvec.clear();
      hets_readc_pvec.clear(); num_covg_hets=0;

      for( unsigned k=0; k<num_rnaseq_samples; k++ )
        if( funct_coverage_test( (*snp_it), k ) )  // If it abides current read depth standard
          {
          hets_pval_pvec.push_back( &(*snp_it).sample_vec[k].pval );
          hets_imba_pvec.push_back( &(*snp_it).sample_vec[k].imba );
          hets_readc_pvec.push_back( &(*snp_it).sample_vec[k].totalcount );
          num_covg_hets++;
          }

      if( num_covg_hets >= min_num_hets )
        {  // Calculate Zval and if consistent or not
        zvalue = funct_calc_zvalue( hets_pval_pvec, hets_imba_pvec, hets_readc_pvec );

        if( zvalue > -100 ) rnd_zval_vec.push_back( zvalue );
        }
      }
  }

// Calculate and writeout the pvals corresponding to 2.5, 5 and 10% FDRs, and the full pval hist
cout << "\nWriting the QQplot\n";
funct_calc_qqplot_and_writeout( rnd_zval_vec, obs_zval_vec, rnd_zval_out, qqplot_zval_out, 2 );

qqplot_zval_out.file.flush(); median_out.file.flush(); rnd_zval_out.file.flush();
cout << "\nDone." << endl;

data_in.close(); qqplot_zval_out.close();
median_out.close(); rnd_zval_out.close();
}





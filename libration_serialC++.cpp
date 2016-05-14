/* Author: Neha Jothi
Description: This is the C++ implementation of the serial Perl Libration program written by Kat Volk 
Input: inputs.txt 
Output: securely-resonant-test-particles
        partially-resonant-test-particles
        gnu plot of resonance, which is currently disabled
*/

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<unordered_map>
#include<list>

using namespace std;

// Names of output files
#define SECURE "securely-resonant-test-particles"
#define INSECURE "partially-resonant-test-particles"
#define pi  3.141592653589793f
#define MAXNSECTIONS 500

ofstream o_s(SECURE, ios::binary);
ofstream o_in(INSECURE, ios::binary);

// Input Parameters
int ntp;
float pmax;
string plfile;
float nsections;
int plots;

//Number of rows in each test particle file
const int n_tpp = 9629;

// Global variables
typedef struct{
    int m, n, p, q, r, s;
    float res_angle[n_tpp];
    float ratios[2500];// list of ratios (p/q) that have been checked
    int ratios_n;
}ivs;

// from particles' file
typedef struct {
    float delta_a;
    float lambda_tp[n_tpp];
    float Omega[n_tpp];
    float omega[n_tpp];
    float M[n_tpp];
    float abar;
    float ebar;
    float ibar;
}p_ps;

// from planet file
float lambda_p[n_tpp];
float Omega_p[n_tpp];
float omega_p[n_tpp];
float M_p[n_tpp];
float time_p[n_tpp];
float apl;

void checkers(p_ps &params, int k);
bool checked(ivs &iVars, float ratio);
void checklib(ivs &iVars, p_ps &params, int &rflag, float &rescent,
              float &resamp, int &rsec_total, int k);

void plot_angle(ivs &iVars);
void tp_read(p_ps &params, int k);

int main()
{
    // reading input file     
    ifstream in("inputs.txt");
    for(int s = 0 ; (in>>ntp) ; s++)
    {
		in>>pmax>>plfile>>nsections>>plots;


		cout<< "ntp: "<<ntp << " pmax : " <<pmax << " plfile: " <<plfile
				<< " nsections: " <<nsections<< " plots: "<<plots<<"\n";
    }

    // reading plfile file (planet_read)    
    {
    float atot = 0;
    float l = 0;
    float data[9];
  
    ifstream in(plfile.c_str(), ios::binary);
    
    int s;
    for( s = 0 ; (in>>data[1]) ; s++) //9629
    {
        in>>data[2]>>data[3]>>data[4]
          >>data[5]>>data[6]>>data[7]>>data[8];

        time_p[s]     = data[2];
        atot         += data[3];
        Omega_p[s]    = data[6];
        omega_p[s]    = data[7];
        M_p[s]        = data[8];

        l = Omega_p[s] + omega_p[s] + M_p[s];
	lambda_p[s] = l;  // lambda_p
    }//End of for
 
    apl = atot / float(s);    
    in.close();
    }

    p_ps params;
    // Loop through particles
    for(int k = 1 ; k <= ntp ; k++)
    {
        // reading one particle  (tp_read)
        tp_read(params, k);
        cout<< " Executing particle " << k << "\n";
        // sub checkres
        if( params.delta_a < 3) //otherwise it won't be long-term resonant
        {
	    checkers(params, k);
        }// end of if checkers
    }// end of loop for particles
    return 0;
}

//
void checkers( p_ps &params, int k)
{
    ivs iVars;
    iVars.ratios_n = 0;
    float ratio;
    string rid("");
    float center=0, amplitude=0;
    float fraction;
    int flag=0;
    //Additional level odf parallelism may be possible here    #pragma omp parallel for
    for(iVars.p = 1; iVars.p <= pmax ; iVars.p++) 
    {
        for(iVars.q = iVars.p ; ( iVars.p-iVars.q < pmax && iVars.q > 0) ; iVars.q--)
        {
            ratio = float(iVars.p)/float(iVars.q);
            if( checked(iVars , ratio ) )
            {
                continue;       // We have checked this ratio, go to next q
            } 
            // if this ratio p/q hasn't already been checked
            iVars.ratios[iVars.ratios_n] = ratio;
            iVars.ratios_n++;
            float a_res = pow(ratio, float(2.0/3.0) ) * apl;
            if( (abs( float(float(params.abar - a_res) / a_res) )) > 0.05f  ) 
            {
                continue;
            }
            //#test particle is within 5% of the resonance location
            //#run through all possible resonance angles for that ratio
            for ( iVars.m = iVars.p-iVars.q ; iVars.m >= 0; iVars.m--) 
            {
                for ( iVars.n = iVars.p-iVars.q-iVars.m ; iVars.n >= 0; iVars.n--)
                {
                    for ( iVars.r = iVars.p-iVars.q-iVars.m-iVars.n ;
                          iVars.r >= 0; iVars.r--) 
                    {
                        iVars.s = iVars.p-iVars.q-iVars.m-iVars.n-iVars.r;
                        //#string that will id the results of this specific resonance angle
                        char t_rid[50];
                        int cc = sprintf (t_rid, "%d:%d:%d:%d:%d:%d",iVars.p,iVars.q,iVars.m,iVars.n,iVars.r,iVars.s );
                        t_rid[cc] = 0;
                        rid = t_rid;
                        int rsec_total = 0;
                        checklib(iVars,params, flag, center, amplitude, rsec_total, k);
                        fraction = float( rsec_total / float(nsections) );
                        if( flag )
                        {
                            o_s<<k<<" "<<params.abar<<" "<<params.ebar
                               <<" "<<params.ibar<<" "
                               <<rid<<" "<<fraction<<" "
                               <<center<<" "<<amplitude<<"\n";
                            if( plots )//Currently disabled
                            {
                                plot_angle( iVars );
                                char command[50];
                                int cc = sprintf (command, 
                                      "mv angle.ps %d-%s.ps", k, rid.c_str() );
                                command[cc] = 0;
                                system(command);
                            }
                            return;
                        }
                        else if( fraction > 0.25)
                        {
                            o_in<<k<<" "<<params.abar<<" "<<params.ebar
                                <<" "<<params.ibar<<" "
                               <<rid<<" "<<fraction<<" "
                               <<center<<" "<<amplitude<<"\n";
                            if( plots )
                            {
                                plot_angle( iVars );
                                char command[50];
                                int cc = sprintf (command, 
                                     "mv angle.ps %d-partial-%s.ps", k, rid.c_str() );
                                command[cc] = 0;
                                system(command);
                            }
                            return;
                        }
                    }// end loop for r 
                }// end loop for n
            }// end loop for m
        }// end loop for q
    }// end loop for p
}

//
void checklib(ivs &iVars, p_ps &params, int &rflag, float &rescent, 
              float &resamp, int &rsec_total, int k)
{
    int points        = n_tpp; // number of time points
    int pts_per_sect  = points / nsections; // points per section
    int outlier_limit = 0; 
    int ints[] = {0,90,180,270};
    list<int> lowbnds (ints,ints+4);
	
    float phi = 0;   
    rflag = 0;
    unordered_map<int,int> n_out;
    unordered_map<int,int> nressec;
    int sec_rflag[MAXNSECTIONS] = {0};
    //section resonance flag: 1 if section is res, 0 if not
    unordered_map<int,float> phi_sum, phi_min, phi_max;
    resamp=0, rescent=0; 
    for (list<int>::iterator lowbnd = lowbnds.begin();
         lowbnd != lowbnds.end(); ++lowbnd)
    {
        nressec[*lowbnd] = 0;       
        phi_min[*lowbnd] = 10000;   
        phi_max[*lowbnd] = -10000;  
        phi_sum[*lowbnd] = 0;      
    }

    //calculate resonance angle for all the timepoints
    for(int t = 0 ; t <  points ; t++)
    {
        phi =  iVars.p*params.lambda_tp[t] - iVars.q*lambda_p[t] - 
               iVars.m*(params.Omega[t] + params.omega[t]) -
               iVars.n*params.Omega[t] - iVars.r*(Omega_p[t] + 
               omega_p[t]) - iVars.s*Omega_p[t];
        phi = float(phi*180/pi);
        while( phi > 360.0) { phi = phi - 360.0f; }
        while( phi < 0.0)   { phi = phi + 360.0f; }
        iVars.res_angle[t] = phi;
    }


    for(int sec = 0 ; sec < nsections ; sec++) 
    {
        for (list<int>::iterator lowbnd = lowbnds.begin();
             lowbnd != lowbnds.end(); ++lowbnd)
            n_out[*lowbnd] = 0;       // reset n_out
		
        for(int j = 0 ; j < pts_per_sect ; ++j)
        {
            int t = sec*pts_per_sect + j;
            phi = iVars.res_angle[t];
            for (list<int>::iterator lowbnd = lowbnds.begin();
                 lowbnd != lowbnds.end(); ++lowbnd)
            {
                if( phi < *lowbnd) { phi += 360;}
                if( phi < *lowbnd+5.0 ||  phi > *lowbnd+355.0) { n_out[*lowbnd]++;}
                phi_sum[*lowbnd] += phi;
                if( phi < phi_min[*lowbnd] ){ phi_min[*lowbnd] = phi;}
                if( phi > phi_max[*lowbnd] ){ phi_max[*lowbnd] = phi;}
            }
        }
        for (list<int>::iterator lowbnd = lowbnds.begin();
               lowbnd != lowbnds.end(); ++lowbnd)
        {
            if(n_out[*lowbnd] <= outlier_limit)
            {
                nressec[*lowbnd]++; 
                if( sec_rflag[ sec ] == 0) //this ensures no double counting
                {
                    sec_rflag[sec] = 1; //#flag this section as resonant
                    rsec_total++; //#increment total # of res sections
                }
            }

            n_out[*lowbnd] = 0; 

  		  } // end of loop for lowbnd
    }// end of loop for sec

    for (list<int>::iterator lowbnd = lowbnds.begin();
               lowbnd != lowbnds.end(); ++lowbnd)
    {
        if( !rflag ) 
        {
            float cent = phi_sum[*lowbnd] / float(points + 19);
            while( cent > 360.0 ) {cent -= 360.0;}
            float amp = ( phi_max[*lowbnd] - phi_min[*lowbnd] ) / 2.0f;
            if( amp < 170.0f ) 
            { 
                resamp  = amp;
                rescent = cent;
                rflag   = 1;
            }
            else if( nressec[*lowbnd] > nsections/2 ) 
            {                                      
                cout<<k<<""<<iVars.p<<":"<<iVars.q<<":"<<iVars.m
                    <<":"<<iVars.n<<":"<<iVars.s<<" for "
                    <<nressec[*lowbnd]<<" sections, "<<"lowbnd="<<*lowbnd<<"\n";
                resamp  = 999; 
                rescent = cent;
                rflag   = 1;
            }
        } //#endof if !rflag
    } //#end of loop for lowbnd

}

void tp_read(p_ps &params, int k)
{
    float atot = 0;
    float etot = 0;
    float itot = 0;
    float l = 0;
    float data[9];
    float min_a = 10000;
    float max_a = 0;

    char t_tpfile[50];
    int cc = sprintf (t_tpfile, "aei.%d", k );
    t_tpfile[cc] = 0;  
    string tpfile(t_tpfile);
    ifstream in( tpfile.c_str() );

    int s;
    for( s = 0 ; in>>data[1] ; s++)
    {

        in>>data[2]>>data[3]>>data[4]>>data[5]>>data[6]>>data[7]>>data[8];
        if( data[3] < min_a )
        {
            min_a = data[3];
        }
        if( data[3] > max_a )
        {
            max_a = data[3];
        }
        atot         += data[3];
        etot         += data[4]; 
        itot         += data[5]; 
        params.Omega[s]    = data[6];
        params.omega[s]    = data[7];
        params.M[s]        = data[8];

        l = params.Omega[s] + params.omega[s] + params.M[s];
	params.lambda_tp[s] = l;  

    }//End of loop for reading one particle's file

    params.delta_a  = float(max_a - min_a);
    params.abar     = float(atot / float(s) );
    params.ebar     = float(etot / float(s) );
    params.ibar     = float(itot / float(s) );
}

bool checked(ivs &iVars, float ratio)
{
    for(int i = 0 ; i < iVars.ratios_n ; ++i )
        if(iVars.ratios[i] == ratio)
            return true;
    return false;
}

void plot_angle(ivs &iVars)
{
    int smax = n_tpp;

    ofstream tp_out("temp.txt", ios::binary);

    for(int s = 0 ; s < smax ; s++)
    {
        tp_out<<time_p[s]<<" "<<iVars.res_angle[s]<<"\n";
    }

    system("gnuplot plot-angle.gnu");
    system("rm temp.txt");
}


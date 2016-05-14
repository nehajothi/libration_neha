/*

Serial without global variables
+.75 hours

*/

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<unordered_map>
#include<list>
#include<mpi.h>

using namespace std;

// Names of output files
#define SECURE "securely-resonant-test-particles"
#define INSECURE "partially-resonant-test-particles"

ofstream o_s(SECURE, ios::binary);
ofstream o_in(INSECURE, ios::binary);

#define pi  3.141592653589793f
#define MAXNSECTIONS 500
// Important Parameters
int ntp;
float pmax;
string plfile;
float nsections;
int plots;

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

int main(int argc, char *argv[])
{
    // reading input file     
    {
    ntp = 30;
    pmax = 30;
    plfile = "neptune_aei";
    nsections = 50;
    plots = 0;
    }//End of section

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
 
    apl = atot / float(s);        //Avereging

//    cout<<"\n\n"<<apl<<"\n"<<atot<<"\n"<<s<<"\n\n";
//    printf("\n%5.9f\n%5.6f\n%d\n\n",apl,atot,s);
    in.close();
    }//End of section
    int myrank,numnodes;
    double startTime;
    MPI_Init(&argc, &argv);
  
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);
    p_ps params;
    // Loop through particles
    if (myrank == 0) {
        startTime = MPI_Wtime();
    }
    for(int k = ((myrank * ntp)/numnodes )+1 ; k <= (myrank + 1)*ntp/numnodes ; k++) //ntp
    {
        cout<<"rank "<< myrank <<"is executing file "<< k <<"\n";
        // reading one particle  (tp_read)
        tp_read(params, k);
    
        // sub checkres
        if( params.delta_a < 3) //otherwise it won't be long-term resonant
        {
	    checkers(params, k);
        }// end of if checkers

    }// end of loop for particles
    MPI_Finalize();
    cout <<" Total time taken "<< MPI_Wtime() - startTime << "\n";
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

    for(iVars.p = 1; iVars.p <= pmax ; iVars.p++) 
    {
        for(iVars.q = iVars.p ; ( iVars.p-iVars.q < pmax && iVars.q > 0) ; iVars.q--)
        {
            ratio = float(iVars.p)/float(iVars.q);
//if(p==3) cout<<"\nhi q = "<<q<<"  ratio = "<<ratio<<"\n";
            if( checked(iVars , ratio ) )
            {
                continue;       // We have checked this ratio, go to next q
            } 
            // if this ratio p/q hasn't already been checked
            iVars.ratios[iVars.ratios_n] = ratio;//push(@ratios, $ratio);
            iVars.ratios_n++;

            float a_res = pow(ratio, float(2.0/3.0) ) * apl;
//cout<<"\nk = "<<k<<"  p = "<<p<<"  q = "<<q<<"  ares = "<<(abs( float(float(abar - a_res) / a_res) ))<<"   nr = "<<ratios_n<<"\n";
            if( (abs( float(float(params.abar - a_res) / a_res) )) > 0.05f  ) 
            {
//cout<<"\nhi"<<(abs( float(float(abar - a_res) / a_res) ))<<"\n";
                continue;
            }
//cout<<"\n2 k = "<<k<<"  p = "<<p<<"  q = "<<q<<"  ares = "<<(abs( float(float(abar - a_res) / a_res) ))<<"\n";
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
                        int cc = sprintf (t_rid, "%d:%d:%d:%d:%d:%d", 
                                                iVars.p,iVars.q,iVars.m,
                                                iVars.n,iVars.r,iVars.s );
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
                            if( plots )
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
    int outlier_limit = 0; // this can be changed to allow for some outlying points
    //accounts for libration around values other than 180
    int ints[] = {0,90,180,270};
    list<int> lowbnds (ints,ints+4);
	
    float phi = 0;   //resonant argument
    rflag = 0; //flag for the object being resonant
    //number of points outside the bounds for a given lowbnd
    unordered_map<int,int> n_out;
    //number resonant sections for this p,q,m,n,r,s and a given lowbnd)
    unordered_map<int,int> nressec;
    int sec_rflag[MAXNSECTIONS] = {0};
    //section resonance flag: 1 if section is res, 0 if not
    //sum, min, max of phi for a given lowbnd
    unordered_map<int,float> phi_sum, phi_min, phi_max;
    resamp=0, rescent=0; //resonance amplitude, libration center

    for (list<int>::iterator lowbnd = lowbnds.begin();
         lowbnd != lowbnds.end(); ++lowbnd)
    {
        nressec[*lowbnd] = 0;       // reset nressec
        phi_min[*lowbnd] = 10000;   // reset min plaeholders
        phi_max[*lowbnd] = -10000;  // reset max plaeholders
        phi_sum[*lowbnd] = 0;       // reset sum plaeholders
    }

    //calculate resonance angle for all the timepoints
    for(int t = 0 ; t <  points ; t++)
    {
        //calculate the resonance angle
        phi =  iVars.p*params.lambda_tp[t] - iVars.q*lambda_p[t] - 
               iVars.m*(params.Omega[t] + params.omega[t]) -
               iVars.n*params.Omega[t] - iVars.r*(Omega_p[t] + 
               omega_p[t]) - iVars.s*Omega_p[t];

        //all the angles were in radians, but we'll do the 
        phi = float(phi*180/pi);
        // check in degrees

        //get phi in range [0,360]
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
            // get phi in range [lowbnd,lowbnd+360]
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
                nressec[*lowbnd]++; //#increment # of res sections for this lowbnd 
                if( sec_rflag[ sec ] == 0) //this ensures no double counting
                {
                    sec_rflag[sec] = 1; //#flag this section as resonant
                    rsec_total++; //#increment total # of res sections
//if($rsec_total{$rid}>$nsections) {warn "rsec_total is > $nsections?!? >:-(\n"} 
//#just in case double counting occurs
                }
            }

            n_out[*lowbnd] = 0; //reset n_out

  		  } // end of loop for lowbnd
    }// end of loop for sec

    for (list<int>::iterator lowbnd = lowbnds.begin();
               lowbnd != lowbnds.end(); ++lowbnd)
    {
        if( !rflag ) //#so we exit as soon as we've determined resonance
        {
            float cent = phi_sum[*lowbnd] / float(points + 19);
            while( cent > 360.0 ) {cent -= 360.0;}
            float amp = ( phi_max[*lowbnd] - phi_min[*lowbnd] ) / 2.0f;
//cout<<"\n\n AMP = "<<amp<<"  low = "<<*lowbnd<<"  M = "<<phi_max[*lowbnd]<<" MI = "<<phi_min[*lowbnd]<<"\n\n";
            if( amp < 170.0f ) //we have well behaved, continuous 
            {                  //libration around a single libration center
                resamp  = amp;
                rescent = cent;
                rflag   = 1;
            }
            else if( nressec[*lowbnd] > nsections/2 ) //#not as clear cut, 
            {                                        //#but still have a 
                //dominant res center with libration half the time
                cout<<k<<""<<iVars.p<<":"<<iVars.q<<":"<<iVars.m
                    <<":"<<iVars.n<<":"<<iVars.s<<" for "
                    <<nressec[*lowbnd]<<" sections, "<<"lowbnd="<<*lowbnd<<"\n";
                resamp  = 999; //#amplitude can't be well defined
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

    char t_tpfile[50];//="aei.1";
    int cc = sprintf (t_tpfile, "aei.%d", k );
    t_tpfile[cc] = 0;  
//cout<<t_tpfile;
    string tpfile(t_tpfile);
    ifstream in( tpfile.c_str() );

    int s;
    for( s = 0 ; in>>data[1] ; s++) //9629
    {

        in>>data[2]>>data[3]>>data[4]
          >>data[5]>>data[6]>>data[7]>>data[8];

//if(s==0) cout<<data[1]<<"  "<<data[2]<<"  "<<data[3]<<"  "<<data[4]<<"  "<<data[5]<<"  "<<data[6]<<"  "<<data[7]<<"  "<<data[8];

        if( data[3] < min_a )
        {
            min_a = data[3];
        }
        if( data[3] > max_a )
        {
            max_a = data[3];
        }
 
        atot         += data[3]; //for a averaging
        etot         += data[4]; //for a averaging
        itot         += data[5]; //for a averaging
        params.Omega[s]    = data[6];
        params.omega[s]    = data[7];
        params.M[s]        = data[8];

        l = params.Omega[s] + params.omega[s] + params.M[s];
	params.lambda_tp[s] = l;  // lambda_p

    }//End of loop for reading one particle's file

    params.delta_a  = float(max_a - min_a);
    params.abar     = float(atot / float(s) );
    params.ebar     = float(etot / float(s) );
    params.ibar     = float(itot / float(s) );

//    cout<<"\n\n"<<apl<<"\n"<<atot<<"\n"<<s<<"\n\n";
//    printf("\n%5.9f\n%5.9f\n%5.9f\n%5.9f\n%5.9f\n\n",atot,abar,ebar,ibar,delta_a);
}

bool checked(ivs &iVars, float ratio)
{
    for(int i = 0 ; i < iVars.ratios_n ; ++i )
        if(iVars.ratios[i] == ratio)
            return true;

//    ratios[ratios_n] = ratio;
//    ratios_n++;

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


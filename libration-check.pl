#!/usr/bin/perl -w
use List::Util qw(first max min sum);
use POSIX;

#inputs
#first line: number of test particle files to read in
#second line: maximum p for p:q resonances to check
#third line: name of the file containing the planet's orbital evolution
#fourth line: number of windows for the resonance check
#fifth line: produce plots? (1 for yes, 0 for no), requires gnuplot

$ntp = <stdin>;
chomp($ntp);

$pmax = <stdin>;
chomp($pmax);

$plfile = <stdin>;
chomp($plfile);

$nsections = <stdin>;
chomp($nsections);

$plots = <stdin>;
chomp($plots);

##################################################################
#
# read in the planet's orbital history
&planet_read($plfile);

open(SECURE,">securely-resonant-test-particles\n");
open(INSECURE,">partially-resonant-test-particles\n");

# loop through the test particle files
tploop: for($k=1;$k<=$ntp;$k++)
	{
	$tpfile = "aei.$k";
	#read in the test particle's orbital history
	&tp_read($tpfile);
	if($delta_a < 3) #otherwise it won't be long-term resonant
		{
		#run through the resonance check
		&checkres;
		}
	}


close(SECURE);
close(INSECURE);




















##################################################################
#subroutines
##################################################################
#
#
#

sub checkres
	{
	@ratios = (); # list of ratios (p/q) that have been checked
	push(@ratios,0);
	for($p=1;$p<=$pmax;$p++) 
		{
	 	for($q=$p;($p-$q<$pmax && $q>0);$q--)
			{
			$ratio = $p/$q;
	 		if(!&checked($ratio)) 
				{ 
				# if this ratio p/q hasn't already been checked
				push(@ratios, $ratio);
				$a_res = $ratio**(2.0/3.0)*$apl;
	 			if( abs(($abar-$a_res)/$a_res) < 0.05 ) 
					{
					#test particle is within 5% of the resonance location
					#run through all possible resonance angles for that ratio
	 				for ($m=$p-$q; $m>=0; $m--) 
						{
	 					for ($n=$p-$q-$m; $n>=0; $n--) 
							{
	 						for ($r=$p-$q-$m-$n; $r>=0; $r--) 
								{
								$s = $p-$q-$m-$n-$r;
								$rid = "$p:$q:$m:$n:$r:$s"; #string that will id the results of this specific resonance angle
 								$rsec_total{$rid} = 0; #number of resonant sections for this resonance angle
								($flag, $center, $amplitude) = &checklib;
								$fraction = $rsec_total{$rid}/$nsections;
								if($flag)
									{
									print SECURE "$k $abar $ebar $ibar $rid $fraction $center $amplitude\n";
									if($plots)
										{
										&plot_angle;
										$plotfile = "$k"."-$rid".".ps";
										`mv angle.ps $plotfile`;
										}
									return;
									}
								elsif($fraction > 0.25)
									{
									print INSECURE "$k $abar $ebar $ibar $rid $fraction $center $amplitude\n";
									if($plots)
										{
										&plot_angle;
										$plotfile = "$k"."-partial-$rid".".ps";
										`mv angle.ps $plotfile`;
										}
									return;

									}
								}
							}
						}
					}
				}
			}
		}
	
	
	 
	 return;
	 
	 
	 
	 
	 
	 }



sub checklib
	{
	my $points = @lambda_tp; #number of time points
 	my $pts_per_sect = $points / $nsections; #points per section
	my $outlier_limit = 0; #this can be changed to allow for some outlying points 
	my $pi = 3.141592653589793;
 	my $sec=0;
   my	$j=0;
	my $t=0; 
 	my @lowbnds=(0,90,180,270); #accounts for libration around values other than 180
	@res_angle = ();
	
 	my $phi=0; #resonant argument
 	my $rflag=0; #flag for the object being resonant
 	my %n_out; #number of points outside the bounds for a given lowbnd
 	my %nressec; #number resonant sections for this p,q,m,n,r,s and a given lowbnd)
 	my @sec_rflag = ();
	for($j=0;$j<$nsections;$j++){push @sec_rflag, 0;}
	#section resonance flag: 1 if section is res, 0 if not
 	my (%phi_sum, %phi_min, %phi_max); #sum, min, max of phi for a given lowbnd
 	my $resamp=0; $rescent=0; #resonance amplitude, libration center

	foreach $lowbnd (@lowbnds) {$nressec{$lowbnd} = 0;} # reset nressec
	foreach $lowbnd (@lowbnds) {$phi_min{$lowbnd} = 1000;} # reset min and max plaeholders
	foreach $lowbnd (@lowbnds) {$phi_max{$lowbnd} = -1000;}

	#calculate resonance angle for all the timepoints
	for($t=0;$t<$points;$t++)
		{
		#calculate the resonance angle
		$phi = $p*$lambda_tp[$t] - $q*$lambda_p[$t] - $m*($Omega[$t]+$omega[$t]) -
	  	$n*$Omega[$t] - $r*($Omega_p[$t]+$omega_p[$t]) - $s*$Omega_p[$t];

		$phi = $phi*180/$pi; #all the angles were in radians, but we'll do the 
		#check in degrees

      # get phi in range [0,360]
      while($phi>360.0) { $phi = $phi-360.0; }
      while($phi<0.0) { $phi = $phi+360.0; }
		
		$res_angle[$t] = $phi;
		}


 	for($sec=0;$sec<$nsections;$sec++) 
		{

		foreach $lowbnd (@lowbnds) {$n_out{$lowbnd} = 0;} # reset n_out
		
  		for($j=0;$j<$pts_per_sect;$j++)
	  		{
			$t = $sec*$pts_per_sect + $j;
			$phi = $res_angle[$t];
			foreach $lowbnd (@lowbnds)
				{
	 			# get phi in range [lowbnd,lowbnd+360]
	 			if($phi<$lowbnd) {$phi += 360;}
	 			if($phi<$lowbnd+5.0 || $phi>$lowbnd+355.0) {$n_out{$lowbnd}++;}
	 			$phi_sum{$lowbnd} +=$phi;
				if($phi<$phi_min{$lowbnd}){$phi_min{$lowbnd}=$phi;}
	 			if($phi>$phi_max{$lowbnd}){$phi_max{$lowbnd}=$phi;}
				}
  			}
  		foreach $lowbnd (@lowbnds)
		 	{
   		if($n_out{$lowbnd}<=$outlier_limit)
				{
				$nressec{$lowbnd}++; #increment # of res sections for this lowbnd
				if($sec_rflag[$sec] == 0) #this ensures no double counting
			  		{
					$sec_rflag[$sec] = 1; #flag this section as resonant
					$rsec_total{$rid}++; #increment total # of res sections
					if($rsec_total{$rid}>$nsections) {warn "rsec_total is > $nsections?!? >:-(\n"} 
					#just in case double counting occurs
					}
   			}
   		$n_out{$lowbnd} = 0; # reset n_out
  			} #endfoeach lowbnd

		} # endfor sec

 	foreach $lowbnd (@lowbnds)
		{
  		if(!$rflag) #so we exit as soon as we've determined resonance
	  		{
			$cent = $phi_sum{$lowbnd} / ($points+19);
			while($cent>360.0) {$cent -= 360.0;}
			$amp = ( $phi_max{$lowbnd}-$phi_min{$lowbnd} ) / 2.;
			if($amp<170.0) #we have well behaved, continuous libration around a single libration center
				{
				$resamp=$amp;
				$rescent=$cent;
				$rflag=1;
				}
			elsif($nressec{$lowbnd}>$nsections/2) #not as clear cut, 
			#but still have a dominant res center with libration half the time
				{
				printf "$k $p:$q:$m:$n:$r:$s for %i sections, lowbnd=%i\n",$nressec{$lowbnd}, $lowbnd;
				$resamp=999; #amplitude can't be well defined
				$rescent=$cent;
				$rflag=1;
				}
  			} #endif !rflag
 		} #endforeach lowbnd

 	return ($rflag, $rescent, $resamp);
	} #endsub checkres






#
#
sub planet_read
	{
	my $file = $_[0];
  	my ($atot, $l, $s);
	open(DAT,"<$file");
	@lambda_p = ();
	@Omega_p = ();
	@omega_p = ();
	@M_p = ();
	@time = ();
	$atot = 0;
	$s = 0;
	while(<DAT>)
		{ 
		$stuff = $_;
		chomp($stuff);
		@data = split(/\s+/,$stuff);
		$atot+=$data[3]; #for a averaging
		$time[$s] = $data[2];
		$Omega_p[$s] = $data[6];
		$omega_p[$s] = $data[7];
		$M_p[$s] = $data[8];
		$l = $Omega_p[$s]+$omega_p[$s]+$M_p[$s];
		$lambda_p[$s] = $l;
		$s+=1;
		} 
		$apl = $atot/($s);
	return;
	}


sub tp_read
	{
	my $file = $_[0];
  	my ($atot, $etot, $itot, $l, $s);
	open(DAT,"<$file");
	@lambda_tp = ();
	@Omega = ();
	@omega = ();
	@M = ();
	my $min_a = 10000;
	my $max_a = 0;

	$atot = 0;
	$etot = 0;
	$itot = 0;
	$s = 0;
	while(<DAT>)
		{ 
		$stuff = $_;
		chomp($stuff);
		@data = split(/\s+/,$stuff);
		if($data[3]<$min_a){$min_a = $data[3];}
		if($data[3]>$max_a){$max_a = $data[3];}
		$atot+=$data[3]; #for a averaging
		$etot+=$data[4]; #for a averaging
		$itot+=$data[5]; #for a averaging
		$Omega[$s] = $data[6];
		$omega[$s] = $data[7];
		$M[$s] = $data[8];
		$l = $Omega[$s]+$omega[$s]+$M[$s];
		$lambda_tp[$s] = $l;
		$s+=1;
		} 
		$abar = $atot/($s);
		$delta_a = $max_a - $min_a;
		$ebar = $etot/($s);
		$ibar = $itot/($s);

	return;
	}



sub checked 
	{
	# checks to see if the ratio is in @ratios
	# (in which case we don't need to check this resonance)
	local $r;
	foreach $r (@ratios) 
		{
  		if ($r==$_[0]) { return 1; }
 		}
 	return 0;
	}

sub plot_angle
	{
	$smax = @time;
	open(OUT,">temp.txt");
	for($s=0;$s<$smax;$s++)
		{
		print OUT "$time[$s] $res_angle[$s]\n";
		}
	close(OUT);
	`gnuplot plot-angle.gnu`;
	`rm temp.txt`;
	return;
	}


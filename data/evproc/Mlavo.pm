#   Copyright (c) 2007 Boulder Real Time Technologies, Inc.           
#                                                                     
#   This software module is wholly owned by Boulder Real Time         
#   Technologies, Inc. This software may be used freely in any 
#   way as long as the copyright statement above is not removed.

package Mlrichter ;

use lib "$ENV{ANTELOPE}/data/evproc" ;

our @ISA= ( "Magnitude" ) ;

use evproc ;

use strict ;
use warnings ;

use lib "$ENV{ANTELOPE}/data/perl" ;

use Datascope ; 

# following are the Richter correction values as a function of distance

our @mltab = (
	(0.0,            1.4),
	(5.0,            1.4),
	(10.0,           1.5),
	(15.0,           1.6),
	(20.0,           1.7),
	(25.0,           1.9),
	(30.0,           2.1),
	(35.0,           2.3),
	(40.0,           2.4),
	(45.0,           2.5),
	(50.0,           2.6),
	(55.0,           2.7),
	(60.0,           2.8),
	(65.0,           2.8),
	(70.0,           2.8),
	(80.0,           2.9),
	(85.0,           2.9),
	(90.0,           3.0),
	(95.0,           3.0),
	(100.0,          3.0),
	(110.0,          3.1),
	(120.0,          3.1),
	(130.0,          3.2),
	(140.0,          3.2),
	(150.0,          3.3),
	(160.0,          3.3),
	(170.0,          3.4),
	(180.0,          3.4),
	(190.0,          3.5),
	(200.0,          3.5),
	(210.0,          3.6),
	(220.0,          3.65),
	(230.0,          3.7),
	(240.0,          3.7),
	(250.0,          3.8),
	(260.0,          3.8),
	(270.0,          3.9),
	(280.0,          3.9),
	(290.0,          4.0),
	(300.0,          4.0),
	(310.0,          4.1),
	(320.0,          4.1),
	(330.0,          4.2),
	(340.0,          4.2),
	(350.0,          4.3),
	(360.0,          4.3),
	(370.0,          4.3),
	(380.0,          4.4),
	(390.0,          4.4),
	(400.0,          4.5),
	(410.0,          4.5),
	(420.0,          4.5),
	(430.0,          4.6),
	(440.0,          4.6),
	(450.0,          4.6),
	(460.0,          4.6),
	(470.0,          4.7),
	(480.0,          4.7),
	(490.0,          4.7),
	(500.0,          4.7),
	(510.0,          4.8),
	(520.0,          4.8),
	(530.0,          4.8),
	(540.0,          4.8),
	(550.0,          4.8),
	(560.0,          4.9),
	(570.0,          4.9),
	(580.0,          4.9),
	(590.0,          4.9),
	(600.0,          4.9) ) ;

sub compml {
	my $self = shift ;
	my $sta = shift ;
	my $millimeters = shift ;

	my $distance = $self->{stations}{$sta}{delta}*111.11 ;

	if ($distance < 0.0 || $distance > 600.0) {return;}
	if ( $millimeters <= 0.0 ) {return;}
	my $i;
	for ($i=0; $i<scalar(@mltab); $i+=2) {
		if ($distance <= $mltab[$i]) {last;}
	}
	my $ml = log($millimeters)/log(10) + $mltab[$i+1];
	return $ml ;
}

sub new {
	return Magnitude::new @_ ;
}

sub getwftimes {
	my $self = shift ;

	my $ret = setup_processes $self ;

	if ($ret ne "ok" ) { return makereturn ( $self, $ret ) ; }

	$self->{stations} = {} ;

	my ($otime,$odepth,$oauth) = dbgetv ( @{$self->{dbo}}, "time", "depth", "auth" ) ;
	my $date = yearday ( $otime ) ;

	if ( defined $self->{params}{auth_accept} ) {
		my $ok = dbex_eval ( @{$self->{dbo}}, "auth =~ /$self->{params}{auth_accept}/" ) ;
		if ( ! $ok ) {
			addlog ( $self, 1, "wrong origin auth " . $oauth ) ;
			return makereturn ( $self, "skip" ) ; 
		}
	}

	my $event_tend = -1.e20 ;
	for ($self->{dba}[3] = 0; $self->{dba}[3] < $self->{nassoc}; $self->{dba}[3]++) {
		my ($sta, $delta) = dbgetv ( @{$self->{dba}} , "sta", "delta" ) ;
		if ( defined $self->{stations}{$sta} ) { next ; }

		my $process ;
		my $channels = {};
		my $ndbv ;

		($ret, $process, $channels, $ndbv) = match_sta ($self, $sta, $otime) ;
		if ( $ret ne "ok" ) { next; }

		if ($delta*111.1 > 600.0) {
			addlog ( $self, 1, $sta . ": station too far away" ) ;
			next ;
		}

		my $pt = dbex_eval ( @{$self->{dbo}}, "ptime(" . $delta . "," . $odepth . ")" ) ;
		my $st = dbex_eval ( @{$self->{dbo}}, "stime(" . $delta . "," . $odepth . ")" ) ;

		my $twin = $process->{signal_twin} ;
		if ( substr($process->{signal_twin}, 0, 1) eq "f") {
			my $fac = substr($process->{signal_twin}, 1) ;
			$twin = 1.1 * $fac * ($st - $pt) ;
		}

		my $noise_twin = $process->{noise_twin};
		if ($process->{noise_twin} eq "tproc") {
			$noise_twin = $twin ;
			if ($noise_twin > 60.0) {$noise_twin = 60.0 ;}
		}

		my $noise_tstart = $otime + $pt - $noise_twin - $process->{noise_toffset} ;
		my $noise_tend = $noise_tstart + $noise_twin ;
		my $signal_tstart = $otime + $pt - $process->{signal_toffset} ;
		my $signal_tend = $signal_tstart + $twin + 10.0 ;

		my $tstart = $noise_tstart - 100.0 ;
		my $tend = $signal_tend ;

		my $hash = {
			"chan_expr" => $process->{chan_expr},
			"delta" => $delta,
			"tstart" => $tstart,
			"tend"	=> $tend,
			"noise_tstart" => $noise_tstart,
			"noise_tend"	=> $noise_tend,
			"signal_tstart" => $signal_tstart,
			"signal_tend"	=> $signal_tend,
			"noise_twin" => $noise_twin,
			"snr_thresh" => $process->{snr_thresh},
			"tupdate" => $self->{params}{update_time},
			"nchans" => $ndbv,
			"channels" => $channels,
			"disposition" => "DataNotReady",
		} ;
		if ( defined $process->{clip_upper} && defined $process->{clip_lower} ) {
			$hash->{clip_upper} = $process->{clip_upper} ;
			$hash->{clip_lower} = $process->{clip_lower} ;
		}
		if ( defined $process->{filter} ) {
			$hash->{filter} = $process->{filter} ;
			if ($hash->{filter} eq "auto") {
				my $expr = sprintf 
					'sta == "%s" && chan =~ /%s/ && %.3f >= time && ( %.3f <= endtime || endtime == null("endtime") )',
						$sta, $process->{chan_expr}, $otime, $otime ;
				my @dbv = dbsubset ( @{$self->{dbc}}, $expr ) ;
				my $ndbv = dbquery ( @dbv, "dbRECORD_COUNT" ) ;
		
				if ($ndbv < 1) {
					addlog ( $self, 0, "station ". $sta . ": no channel matches to "
								. $process->{chan_expr} . " in calibration table" ) ;
					undef $hash->{filter} ;
				} else {
					$dbv[3] = 0;
					my $segtype = dbgetv (@dbv, "segtype");
					if ($segtype eq "V") {
						$hash->{filter} = "WAV" ;
					} elsif ($segtype eq "A") {
						$hash->{filter} = "WAA" ;
					} else {
						addlog ( $self, 0, "station ". $sta . 
							" Cannot determine auto filter for segtype " . $segtype ) ;
						undef $hash->{filter} ;
					}
				}
				dbfree @dbv ;
			} elsif ($hash->{filter} eq "autosp") {
				my $expr = sprintf 
					'sta == "%s" && chan =~ /%s/ && %.3f >= time && ( %.3f <= endtime || endtime == null("endtime") )',
						$sta, $process->{chan_expr}, $otime, $otime ;
				my @dbv = dbsubset ( @{$self->{dbc}}, $expr ) ;
				my $ndbv = dbquery ( @dbv, "dbRECORD_COUNT" ) ;
		
				if ($ndbv < 1) {
					addlog ( $self, 0, "station ". $sta . ": no channel matches to "
								. $process->{chan_expr} . " in calibration table" ) ;
					undef $hash->{filter} ;
				} else {
					$dbv[3] = 0;
					my $segtype = dbgetv (@dbv, "segtype");
					if ($segtype eq "V") {
						$hash->{filter} = 'INT s0.2;G 2080.0 1.e-6' ;
					} elsif ($segtype eq "A") {
						$hash->{filter} = 'INT2 s0.2;G 2080.0 1.e-6' ;
					} else {
						addlog ( $self, 0, "station ". $sta . 
							" Cannot determine auto filter for segtype " . $segtype ) ;
						undef $hash->{filter} ;
					}
				}
				dbfree @dbv ;
			}
		}
		$self->{stations}{$sta} = $hash ;
		if ( $signal_tend > $event_tend ) { $event_tend = $signal_tend; }
	}

#	display $self ;

	if ( scalar ( keys ( %{$self->{stations}} ) ) < 1 ) {
		addlog ( $self, 0, "No channels to process" ) ;
		return makereturn ( $self, "ok" ) ; 
	}

	if ( defined $self->{params}{maximum_bad_fraction} ) {
		$self->{maximum_bad_fraction} = $self->{params}{maximum_bad_fraction} ;
	} else {
		$self->{maximum_bad_fraction} = 0.0;
	}

	if ( defined $self->{params}{maximum_wait_time} ) {
		$self->{expire_time} = $event_tend + $self->{params}{maximum_wait_time} ;
		my $now_time = now() + $self->{params}{maximum_wait_time} ;
		if ( $now_time > $self->{expire_time} ) {
			$self->{expire_time} = $now_time ;
		}
	}

	if ( defined $self->{expire_time} ) {
		return makereturn ( $self, "ok", "stations" => $self->{stations},
				"expire_time" => $self->{expire_time} ) ;
	} else {
		return makereturn ( $self, "ok", "stations" => $self->{stations} ) ;
	}
}

sub process_channel {
	my $self = shift ;
	my $ret = $self->SUPER::process_channel(@_) ;

	if ( $ret->{disposition} ne "channeldone" 
		&& $ret->{disposition} ne "stationdone"
		&& $ret->{disposition} ne "processdone" ) {return $ret;}

	my $sta = $ret->{sta} ;
	my $chan = $ret->{chan} ;

	if ( defined $self->{stations}{$sta}{channels}{$chan}{is_nullcalib} 
			&& $self->{stations}{$sta}{channels}{$chan}{is_nullcalib} ) {
		addlog ( $self, 1, "%s: %s: Channel mag not computed because of null calib",
 						$sta, $chan )  ;
		$self->{stations}{$sta}{disposition} = "NullCalib" ;
		return $ret ;
	}
	if ( defined $self->{stations}{$sta}{channels}{$chan}{is_clipped} 
			&& $self->{stations}{$sta}{channels}{$chan}{is_clipped} ) {
		addlog ( $self, 1, "%s: %s: Channel mag not computed because of clipped data",
 						$sta, $chan )  ;
		$self->{stations}{$sta}{disposition} = "DataClipped" ;
		return $ret ;
	}
	if ( ! defined $self->{stations}{$sta}{channels}{$chan}{signal_amax} ) {
		addlog ( $self, 1, "%s: %s: Channel mag not computed because of no data",
 						$sta, $chan )  ;
		return $ret ;
	}
 	if ( defined $self->{stations}{$sta}{channels}{$chan}{snr} ) {
		if ( $self->{stations}{$sta}{snr_thresh} < 1.0
				|| $self->{stations}{$sta}{channels}{$chan}{snr}
					> $self->{stations}{$sta}{snr_thresh} ) {
			my $millimeters =
 				$self->{stations}{$sta}{channels}{$chan}{signal_amax} ;
			if ( $self->{stations}{$sta}{snr_thresh} >= 1.0 ) {
				$millimeters -= 
 					$self->{stations}{$sta}{channels}{$chan}{noise_std} ;
			} 
 			$self->{stations}{$sta}{channels}{$chan}{m} = compml ( 
				$self, $sta, $millimeters ) ;
 			$self->{stations}{$sta}{channels}{$chan}{m_time} = 
				$self->{stations}{$sta}{channels}{$chan}{signal_tmax} ;
 			$self->{stations}{$sta}{channels}{$chan}{m_snr} = 
				$self->{stations}{$sta}{channels}{$chan}{snr} ;
 			$self->{stations}{$sta}{channels}{$chan}{m_val1} = $millimeters ;
 			$self->{stations}{$sta}{channels}{$chan}{m_units1} = "mmwa" ;
			addlog ( $self, 1, "%s: %s: Channel mag = %.3f",
 					$sta, $chan,
 					$self->{stations}{$sta}{channels}{$chan}{m} ) ;
		} else {
			addlog ( $self, 1, "%s: %s: Channel mag not computed because of low snr",
 							$sta, $chan )  ;
			$self->{stations}{$sta}{disposition} = "LowSnr" ;
				
		}
	} else {
 		$self->{stations}{$sta}{channels}{$chan}{m} = compml ( $self, $sta, 
 				$self->{stations}{$sta}{channels}{$chan}{signal_amax} ) ;
 		$self->{stations}{$sta}{channels}{$chan}{m_time} = 
				$self->{stations}{$sta}{channels}{$chan}{signal_tmax} ;
 		$self->{stations}{$sta}{channels}{$chan}{m_snr} = 
				$self->{stations}{$sta}{channels}{$chan}{snr} ;
 		$self->{stations}{$sta}{channels}{$chan}{m_val1} = 
 				$self->{stations}{$sta}{channels}{$chan}{signal_amax} ;
 		$self->{stations}{$sta}{channels}{$chan}{m_units1} = "mmwa" ;
		addlog ( $self, 2, "%s: %s: Channel mag = %.3f",
 					$sta, $chan,
 					$self->{stations}{$sta}{channels}{$chan}{m} ) ;
	}

	return $ret ;
}

sub process_network {
	my $self = shift ;
	#my $ret = $self->SUPER::process_network(@_) ;
	my $ret = $self->process_network_MAGNITUDE(@_) ;

	my @dborigin = @{$self->{dbo}} ;
	$dborigin[3] = 0 ;
	my $auth = dbgetv ( @dborigin, "auth" ) ;
	dbputv ( @dborigin, "auth", $auth . "Ml" ) ;

	if (defined $self->{m_median} ) {
		dbputv ( @dborigin, "ml", $self->{m_median}, "mlid", $self->{magid} ) ;
	}

	return $ret ;
}

sub process_network_MAGNITUDE {
        my $self = shift ;
        my $flush = shift ;

        my $disp = "ok" ;

        undef $self->{output}{db} ;

        my $m = 0.0;
        my $nm = 0;
        my $std = 0;

        my @mags ;
        my $sta ;
        my @logs ;
        my $log ;
        my $nstas = 0 ;
        foreach $sta (keys(%{$self->{stations}})) {
                $nstas ++ ;
                if ( defined $log ) {
                        $log .= ', ' ;
                }
                $log .= $sta . "=" . $self->{stations}{$sta}{disposition}  ;
                if ( length $log > 80 ) {
                        push @logs, $log ;
                        undef $log ;
                }
                if ( ! defined $self->{stations}{$sta}{m} ) { next; }
                $m += $self->{stations}{$sta}{m} ;
                $std += $self->{stations}{$sta}{m} * $self->{stations}{$sta}{m} ;
                $nm ++ ;
                push @mags, $self->{stations}{$sta}{m} ;
        }

        if ( defined $log ) {
                push @ logs, $log ;
        }

        if ( @logs ) {
                foreach $log ( @logs ) {
                        addlog ( $self, 1, "Network mag: " . $log ) ;
                }
        }

        if ( $nm < 1 ) {
                addlog ( $self, 1, "Network mag: No data" ) ;
                return makereturn ( $self, "nodata" ) ;
        }

        $m /= $nm ;
        $std /= $nm ;
        $std -= $m*$m ;
        $self->{m_mean} = $m ;
        $self->{m_n} = $nm ;
        $self->{m_std} = sqrt ( $std ) ;
        $self->{m_pcnt} = 100.0 * $nm / $nstas ;

        my @magss = sort { $a <=> $b } @mags ;
        my $n = scalar ( @magss ) ;
        if ( $n % 2 ) {
                $self->{m_median} = $magss[int($n/2)] ;
        } else {
                $self->{m_median} = 0.5*$magss[int($n/2-1)] ;
                $self->{m_median} += 0.5*$magss[int($n/2)] ;
        }

        my $lo = $magss[int(0.1587*$n)] ;
        my $hi = $magss[int(0.8413*$n)] ;
        $self->{m_unc} = 0.5*($hi-$lo) ;

        addlog ( $self, 1, "Network mag: mean = %.2f, std = %.2f, median = %.2f, unc = +%.2f/-%.2f, n = %d of %d, pcnt = %.1f%%",
                $self->{m_mean}, $self->{m_std}, $self->{m_median},
                $hi-$self->{m_median}, $self->{m_median}-$lo, $self->{m_n}, $nstas, $self->{m_pcnt} ) ;

        if (defined $self->{params}{station_number_minimum}) {
                if ( $self->{m_n} < $self->{params}{station_number_minimum} ) {
                        addlog ( $self, 1, "Network mag: Rejected - Number of stations %d < minimum %d",
                                $self->{m_n},
                                $self->{params}{station_number_minimum} ) ;
                        undef $self->{m_median} ;
                        return makereturn ( $self, "nodata" ) ;
                }
        }

        if (defined $self->{params}{station_percentage_minimum}) {
                if ( $self->{m_pcnt} < $self->{params}{station_percentage_minimum} ) {
                        addlog ( $self, 1, "Network mag: Rejected - Percentage of stations %.1f%% < minimum %.1f%%",
                                $self->{m_pcnt},
                                $self->{params}{station_percentage_minimum} ) ;
                        undef $self->{m_median} ;
                        return makereturn ( $self, "nodata" ) ;
                }
        }

        if (defined $self->{params}{uncertainty_maximum}) {
                if ( abs($hi-$self->{m_median}) > $self->{params}{uncertainty_maximum} ||
                     abs($self->{m_median}-$lo) > $self->{params}{uncertainty_maximum} ) {
                        addlog ( $self, 1, "Network mag: Rejected - Uncertainty = +%.2f/-%.2f > maximum %.2f",
                                abs($hi-$self->{m_median}),
                                abs($self->{m_median}-$lo),
                                $self->{params}{uncertainty_maximum} ) ;
                        undef $self->{m_median} ;
                        return makereturn ( $self, "nodata" ) ;
                }
        }

        my @dbnetmag = dblookup ( @{$self->{db}}, 0, "netmag", "dbALL", "dbNULL" ) ;
        dbget ( @dbnetmag, 0 ) ;
        @dbnetmag = dblookup ( @dbnetmag, 0, "netmag", 0, "dbSCRATCH" ) ;
        my $magid = dbnextid ( @dbnetmag, "magid" ) ;
        dbputv ( @dbnetmag, "orid", $self->{orid}, "evid", $self->{evid},
                                        "magid", $magid,
                                        "magtype", $self->{params}{output_magtype},
                                        "nsta", $self->{m_n},
                                        "magnitude", $self->{m_median},
                                        "uncertainty", $self->{m_unc},
                                        "auth", $self->{params}{output_auth} ) ;
        $self->{magid} = $magid ;
        my $rec = dbadd ( @dbnetmag ) ;
        $dbnetmag[3] = $rec;

        $self->{output}{db}{assoc_params}{smart_assoc} = "yes";
        $self->{output}{db}{assoc_params}{magnitude_update} = "yes";
        push @{$self->{output}{db}{tables}}, $self->{dbo} ;
        push @{$self->{output}{db}{tables}}, \@dbnetmag ;

        my ( @dbstamag, @dbarrival, @dbwfmeas ) ;

        if ( isyes $self->{params}{output_stamag} ) {
                @dbstamag = dblookup ( @{$self->{db}}, 0, "stamag", 0, 0 ) ;
                @dbarrival = dblookup ( @{$self->{db}}, 0, "arrival", 0, 0 ) ;
                @dbwfmeas = dblookup ( @{$self->{db}}, 0, "wfmeas", 0, 0 ) ;
                my ( $stamag0, $arrival0, $wfmeas0 ) ;
                my ( $stamag1, $arrival1, $wfmeas1 ) ;
                foreach $sta (keys(%{$self->{stations}})) {
                        if ( ! defined $self->{stations}{$sta}{m} ) {next; }
                        my $arid = -1 ;
                        if ( isyes $self->{params}{output_wfmeas} ) {
                                @dbarrival = dblookup ( @dbarrival, 0, 0, "dbALL", "dbNULL" ) ;
                                dbget ( @dbarrival, 0 ) ;
                                @dbarrival = dblookup ( @dbarrival, 0, 0, 0, "dbSCRATCH" ) ;
                                $arid = dbnextid ( @dbarrival, "arid" ) ;
                                dbputv ( @dbarrival, "sta", $sta, "chan", $self->{stations}{$sta}{m_chan},
                                                        "arid", $arid,
                                                        "time", $self->{stations}{$sta}{m_time},
                                                        "jdate", yearday($self->{stations}{$sta}{m_time}),
                                                        "iphase", $self->{params}{output_magtype},
                                                        "amp", $self->{stations}{$sta}{m_amp},
                                                        "per", $self->{stations}{$sta}{m_per},
                                                        "logat", $self->{stations}{$sta}{m_logat},
                                                        "snr", $self->{stations}{$sta}{m_snr},
                                                        "auth", $self->{params}{output_auth} ) ;
                                $rec = dbadd ( @dbarrival ) ;
                                $dbarrival[3] = $rec;
                                if ( ! defined $arrival0 ) { $arrival0 = $rec; }
                                $arrival1 = $rec+1;
                                @dbwfmeas = dblookup ( @dbwfmeas, 0, 0, "dbALL", "dbNULL" ) ;
                                dbget ( @dbwfmeas, 0 ) ;
                                @dbwfmeas = dblookup ( @dbwfmeas, 0, 0, 0, "dbSCRATCH" ) ;
                                dbputv ( @dbwfmeas, "sta", $sta, "chan", $self->{stations}{$sta}{m_chan},
                                                "arid", $arid,
                                                "meastype", $self->{params}{output_magtype},
                                                "time", $self->{stations}{$sta}{signal_tstart},
                                                "endtime", $self->{stations}{$sta}{signal_tend},
                                                "tmeas", $self->{stations}{$sta}{m_time},
                                                "twin", $self->{stations}{$sta}{m_twin},
                                                "val1", $self->{stations}{$sta}{m_val1},
                                                "val2", $self->{stations}{$sta}{m_val2},
                                                "units1", $self->{stations}{$sta}{m_units1},
                                                "units2", $self->{stations}{$sta}{m_units2},
                                                "auth", $self->{params}{output_auth} ) ;
                                if (defined $self->{stations}{$sta}{filter}) {
                                        dbputv ( @dbwfmeas,
                                                "filter", $self->{stations}{$sta}{filter} ) ;
                                }
                                $rec = dbadd ( @dbwfmeas ) ;
                                if ( ! defined $wfmeas0 ) { $wfmeas0 = $rec; }
                                $wfmeas1 = $rec+1;
                        }
                        @dbstamag = dblookup ( @dbstamag, 0, 0, "dbALL", "dbNULL" ) ;
                        dbget ( @dbstamag, 0 ) ;
                        @dbstamag = dblookup ( @dbstamag, 0, 0, 0, "dbSCRATCH" ) ;
                        dbputv ( @dbstamag, "magid", $magid, "sta", $sta, "orid", $self->{orid},
                                                "evid", $self->{evid},
                                                "arid", $arid,
                                                "phase", $self->{params}{output_magtype},
                                                "magtype", $self->{params}{output_magtype},
                                                "magnitude", $self->{stations}{$sta}{m},
                                                "auth", $self->{params}{output_auth} ) ;
                        $rec = dbadd ( @dbstamag ) ;

                        if ( ! defined $stamag0 ) { $stamag0 = $rec; }
                        $stamag1 = $rec+1;

                }
                if ( defined $stamag0 ) {
                        $dbstamag[3] = $stamag0;
                        $dbstamag[2] = $stamag1;
                        push @{$self->{output}{db}{tables}}, \@dbstamag ;
                }
                if ( defined $arrival0 ) {
                        $dbarrival[3] = $arrival0;
                        $dbarrival[2] = $arrival1;
                        push @{$self->{output}{db}{tables}}, \@dbarrival ;
                }
                if ( defined $wfmeas0 ) {
                        $dbwfmeas[3] = $wfmeas0;
                        $dbwfmeas[2] = $wfmeas1;
                        push @{$self->{output}{db}{tables}}, \@dbwfmeas ;
                }

        }

        return makereturn ( $self, $disp ) ;
}




1;

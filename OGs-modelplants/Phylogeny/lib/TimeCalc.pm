package TimeCalc;
use DateTime;

# Retourne l'heure actuelle en secondes
sub temps
{
    # r�cup�ration de l'heure
    my $heure = DateTime->now();
    # Mise au format hh:mm:ss
    my $format = $heure->hms();
    # Split
    my @temp = split(/:/,$format);
    my $secondes = $temp[0]*60*60 + $temp[1]*60 +$temp[2];
    return $secondes;
}

# Retourne le nombre de secondes du temps qui s�pare les deux heures:
# Les entr�es sont en secondes
sub secondes
{
    my $class = shift @_;
    my ($t1,$t2) = @_;
    ($t2-$t1)%60;
}

# Retourne le nombre de minutes du temps qui s�pare les deux heures:
sub minutes
{
    my $class = shift @_;
    my ($t1,$t2) = @_;
    my $sec = ($t2-$t1)%60;
    ($t2-$t1-$sec)/60;
}

1;

# cecile pereira (revu par CD)
# convertir les 2 fichiers de sortie de Inparanoid au format d'entree de MaRiO
# 1 groupe par ligne
# les IDs des proteines du groupe sont séparées par des ;
# les inparalogues sont ajoutés dans les groupes corespondants
use tab2mario;
if ($#ARGV < 0)
{
	print "perl inparanoid2mario.pl <git directory>\n";
	exit();
}

my $work_dir = $ARGV[O];	

$inparanoid_output_dir = "$work_dir/Inparanoid/STEP3/";
$mario_convert_dir = "$work_dir/mario_converted/";
mkdir ($mario_convert_dir) unless (-e $mario_convert_dir);
chdir($inparanoid_output_dir);
tab2mario::convert("inparanoid",$inparanoid_output_dir, $mario_convert_dir);


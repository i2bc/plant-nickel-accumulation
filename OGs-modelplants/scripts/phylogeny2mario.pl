# cecile pereira (revu par CD)
# convertir les 2 fichiers de sortie de Phylogenie au format d'entree de MaRiO
# 1 groupe par ligne
# les IDs des proteines du groupe sont séparées par des ;
# les inparalogues sont ajoutés dans les groupes corespondants
use tab2mario;
if ($#ARGV < 0)
{
	print "perl phylogeny2mario.pl <git directory> [run tag]\n";
	exit();
}

my $work_dir = $ARGV[O];	
my $tag = $ARGV[1];
if($tag) { $tag = '_'.$tag;} else { $tag = ''};

$phylogeny_output_dir = "$work_dir/Phylogeny/STEP5$tag/";
$mario_convert_dir = "$work_dir/mario_converted$tag/";
mkdir ($mario_convert_dir) unless (-e $mario_convert_dir);
chdir($phylogeny_output_dir);
tab2mario::convert("phylogeny",$phylogeny_output_dir, $mario_convert_dir);


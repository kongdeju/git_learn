#!/usr/bin/perl -w
use Getopt::Long;
use Bio::SearchIO;
use File::Basename;
my $q_gz1 = shift;
my $q_gz2 = shift;
my $q_file1=fileparse($q_gz1,'.gz');
my $q_file2=fileparse($q_gz2,'.gz');
my $trim_file1="$q_file1" . "\.trimmed";
my $trim_file2="$q_file2" . "\.trimmed";
my $paired_file1="$trim_file1" . "\.paired1";
my $paired_file2="$trim_file1" . "\.paired2";
my $out_dir=substr($q_gz1,0,7);
##print "$q_gz1 $q_gz2 $q_file1 $q_file2 $trim_file1 $trim_file2 $paired_file1 $paired_file2 $out_dir";
my $s=41;
my $e=91;
my $x=4;
my $t=5;
my $velvetg="-cov_cutoff 20";
my $velveth="-fastq -shortPaired temp.fastq";
GetOptions(	   "s:i"=>\$s,
                   "e:i"=>\$e,
                   "x:i"=>\$x,
                   "t:i"=>\$t,
                   "g:s"=>\$velvetg,
                   "h:s"=>\$velveth,
                   "o:s"=>\$out_dir,
                   );
$Usage="        Usage:  perl Assmbly_Kong.pl [option] <fastq1.gz> <fastq2.gz>\n
        Option:
                s=[int]         the minimum of kerm (default 51);
                e=[int]         the maximum of kerm (default 99);
                x=[int]         the step of minimum kerm to maximum kerm (default 2);
                t=[int]         the max thread (default 10);
                g=[string]      the command option of velvetg (default '-fastq -shortPaired temp.fastq');
                h=[string]      the command option of velveth (default '-cov_cutoff 20');
                o=[string]      the output direction name (default Velvet_out_dir);
                \n";
unless($q_gz1 and $q_gz2){
        print $Usage;
        exit;
}
system "gzip -d $q_gz1 $q_gz2";
system "DynamicTrim.pl $q_file1";
system "DynamicTrim.pl $q_file2";
system "LengthSort.pl $trim_file1 $trim_file2";
system "shuffleSequences_fastq.pl $paired_file1 $paired_file2 temp.fastq";
system "VelvetOpimiser -s $s -e $e -x $x -d $out_dir -t $t -f '$velveth' -o '$velvetg'";
chdir "./$out_dir";
system "perl /work/kondeju/TOOLS/Predict.pl contigs.fa $out_dir -G &";
my $name="$out_dir" . "\.faa";
my $O_db="/work/kongdeju/O_K_database/O_.faa";
`formatdb -i $O_db `;
`blastall -p blastp -i $name -d $O_db -o ${out_dir}_O.bsp -F F -e 1e-5 -b 1 -v 1 -a 10`;
open OUT ,">${out_dir}_O.xls";
my $bsp = Bio::SearchIO->new(-format=>'blast',-file=>"{$out_dir}_O.bsp");
while(my $result = $bsp->next_result){
        my $query_name = $result->query_name;
        my $query_desc = $result->query_description;
        if(my $hit= $result->next_hit){
                my $hit_name = $hit->name;
                my $hit_desc = $hit->description;
                my $hsp = $hit->next_hsp;
                my $ident=$hsp->frac_identical('total')*100;
                my $conse=$hsp->frac_conserved('total')*100;
                print OUT "$query_name\t$hit_name\t$hit_desc\t$ident\t$conse\n";
        }   
        else{
                print OUT "$query_name\t$query_desc\t-\t-\t-\t-\n";
        }   
}
close (OUT);
my $K_db="/work/kongdeju/O_K_database/K_.faa";
`formatdb -i $K_db`;
`blastall -p blastp -i $name -d $K_db -o ${out_dir}_K.bsp -F F -e 1e-5 -b 1 -v 1 -a 10`;
open OUT ,">${out_dir}_K.xls";
my $bsp1 = Bio::SearchIO->new(-format=>'blast',-file=>"${out_dir}_K.bsp");
while(my $result = $bsp1->next_result){
        my $query_name = $result->query_name;
        my $query_desc = $result->query_description;
        if(my $hit= $result->next_hit){
                my $hit_name = $hit->name;
                my $hit_desc = $hit->description;
                my $hsp = $hit->next_hsp;
                my $ident=$hsp->frac_identical('total')*100;
                my $conse=$hsp->frac_conserved('total')*100;
                print OUT "$query_name\t$hit_name\t$hit_desc\t$ident\t$conse\n";
        }   
        else{
                print OUT "$query_name\t$query_desc\t-\t-\t-\t-\n";
        }   
}
close (OUT);
my $T_db="/work/kongdeju/O_K_database/t_.faa";
`formatdb -i $T_db `;
`blastall -p blastp -i $name -d $T_db -o ${out_dir}_T.bsp -F F -e 1e-5 -b 1 -v 1 -a 10`;
open OUT ,">${out_dir}_T.xls";
my $bsp2 = Bio::SearchIO->new(-format=>'blast',-file=>"${out_dir}_T.bsp");
while(my $result = $bsp->next_result){
        my $query_name = $result->query_name;
        my $query_desc = $result->query_description;
        if(my $hit= $result->next_hit){
                my $hit_name = $hit->name;
                my $hit_desc = $hit->description;
                my $hsp = $hit->next_hsp;
                my $ident=$hsp->frac_identical('total')*100;
                my $conse=$hsp->frac_conserved('total')*100;
                print OUT "$query_name\t$hit_name\t$hit_desc\t$ident\t$conse\n";
        }   
        else{
                print OUT "$query_name\t$query_desc\t-\t-\t-\t-\n";
        }   
}
close (OUT);

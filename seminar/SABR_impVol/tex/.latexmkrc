$pdf_mode  = 3;
$dvipdf    = 'dvipdfmx -I 6 %O -o %D %S';
$latex     = 'internal mylatex uplatex  %A %O %S';
$pdflatex  = 'internal mylatex pdflatex %A %O %S';
$bibtex = 'pbibtex %O %B';
$clean_ext = "$clean_ext fmt";

$max_repeat = 5;

sub mylatex {
    my ($engine, $base, @args) = @_;
    my $com = join(' ', @args);
    unless (-e "$base.fmt"){
        print "mylatex: making $base.fmt in ini mode... \n";
        Run_subst("$engine -ini -jobname=\"$base\" \\\&$engine mylatexformat.ltx %S");
    }
    print "mylatex: $base.fmt detected, so running normal latex... \n";
    return Run_subst("$engine -fmt $base $com");
}
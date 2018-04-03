#!/usr/bin/env perl

@default_files=('IAB_NICAM_kernels');
add_cus_dep('glo', 'gls', 0, 'run_makeglossaries');
add_cus_dep('acn', 'acr', 0, 'run_makeglossaries');

sub run_makeglossaries {
  if ( $silent ) {
    system "makeglossaries -q $_[0]";
  }
  else {
    system "makeglossaries $_[0]";
  };
}
push @generated_exts, 'glo', 'gls', 'glg';
push @generated_exts, 'acn', 'acr', 'alg';
$clean_ext .= ' %R.ist %R.xdy';

# add_cus_dep('glo', 'gls', 0, 'run_makeindex');
# add_cus_dep('acn', 'acr', 0, 'run_makeindex');
# sub run_makeindex {
#    my $source = $$Psource;
#    my $dest = $$Pdest;
#    my $log = $dest.".log";
#    my $cmd = "makeindex %O -s \"$_[0].ist\"  -t \"$log\" -o \"$dest\" \"$source\"";
#    if ($silent) { $cmd =~ s/%O/-q/; }
#    else { $cmd =~ s/%O//; }
#    return system $cmd;
# }

if ($^O eq 'MSWin32') {
  $latex = 'latex %O -no-guess-input-enc -synctex=0 -halt-on-error %S';
  $pdflatex = 'pdflatex %O -synctex=1 %S';
  $biber = 'biber %O --bblencoding=utf8 -u -U --output_safechars %B';
  $bibtex = 'bibtex %O %B';
  $makeindex = 'mendex %O -o %D %S';
  $dvipdf = 'dvipdfmx %O -o %D %S';
  $dvips = 'dvips %O -z -f %S | convbkmk -u > %D';
  $ps2pdf = 'ps2pdf.exe %O %S %D';
  $pdf_mode = 1;
  if (-f 'C:/Program Files/SumatraPDF/SumatraPDF.exe') {
    $pdf_previewer = '"C:/Program Files/SumatraPDF/SumatraPDF.exe" -reuse-instance';
  } elsif (-f 'C:/Program Files (x86)/SumatraPDF/SumatraPDF.exe') {
    $pdf_previewer = '"C:/Program Files (x86)/SumatraPDF/SumatraPDF.exe" -reuse-instance';
  } else {
    $pdf_previewer = 'texworks';
  }
} else {
  $latex = 'latex %O -synctex=1 -halt-on-error %S';
  $pdflatex = 'pdflatex %O -synctex=1 %S';
  $biber = 'biber %O --bblencoding=utf8 -u -U --output_safechars %B';
  $bibtex = 'bibtex %O %B';
  $makeindex = 'mendex %O -o %D %S';
  $dvipdf = 'dvipdfmx %O -o %D %S';
  $dvips = 'dvips %O -z -f %S | convbkmk -u > %D';
  $ps2pdf = 'ps2pdf %O %S %D';
  $pdf_mode = 1;
  if ($^O eq 'darwin') {
    $pvc_view_file_via_temporary = 0;
    $pdf_previewer = 'open';
  } else {
    $pvc_view_file_via_temporary = 0;
    $pdf_previewer = 'xdg-open';
  }
}

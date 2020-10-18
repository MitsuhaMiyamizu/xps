use strict;
use File::Find;

# get path containing "..\VC\bin"
my @path = split(';', $ENV{MSVSPATH});
my $pathvc = "";
my $i;
for $i (0..$#path) {
#print ("path: $path[$i]\n"); ####test
   if ($path[$i] =~ /VC\\bin$/) {
      $pathvc = $path[$i];
      last;
   }#if
}#for
#print ("found pathvc: $pathvc\n"); ####test
#print ("found path: @path\n"); ####test

# input file names
my $file2find   = "cl.exe";
my @directories = ($pathvc);

# output file names
my $file_found = "";
my $name_found = "";
my $dir_found  = "";

# check if VC++ is installed in fixed location
find(\&wanted, @directories);

#print ("file_found: $file_found\n"); ####test
#print ("name_found: $name_found\n"); ####test
#print ("dir_found: $dir_found\n"); ####test

# return result
if (!($file_found =~ /^cl.exe$/)) {
   print("ERR_CLEXE");
   exit;
}#if

# check for file vcvarsall.bat
@path      = split('\\\bin', $dir_found);
$file2find = "vcvarsall.bat";
find(\&wanted, @path);

# return result
if ($file_found =~ /^vcvarsall.bat$/) {
#   print ("$name_found");
   print ("$dir_found");
} else {
   print("ERR_VCVARSALLEXE");
}#if


#--- subroutine wanted ---#
sub wanted {
   if ($_ =~ /$file2find$/) {
      $file_found = $_;
      $name_found = $File::Find::name;
      $dir_found  = $File::Find::dir;
      last;
   }#if
}#sub

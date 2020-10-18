use strict;
use File::Find;

# check if ROOTSYS exists
my $rootsys = $ENV{ROOTSYS};
if (!($rootsys =~ /root$/)) {
   print("ERR_ROOTSYS");
   exit;
}#if

# input file names
my $file2find   = "root.exe";
my @directories = ($rootsys);

# output file names
my $file_found;
my $dir_found;

# check if root is installed in fixed location
find(\&wanted, @directories);

sub wanted {
   if ($_ =~ /$file2find$/) { #need to avoid finding "../man1/root.exe.1"
      $file_found = $_;
      $dir_found  = $File::Find::dir;
      last;
   }#if
}#sub

# return result
if ($file_found =~ /^root.exe$/) {
   print ("$file_found");
} else {
   print("ERR_ROOTEXE");
}#if

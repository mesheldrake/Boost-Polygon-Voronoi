#!/usr/bin/perl -w
use strict;
use Module::Build::WithXSpp;

my $build = Module::Build::WithXSpp->new(
  module_name     => 'Boost::Polygon::Voronoi',
  license         => 'perl',
    'build_requires' => {
    'ExtUtils::Typemaps::Default' => '0.05',
    'ExtUtils::XSpp' => '0.16',
    'Module::Build' => '0.3601',
    'Test::More' => '0',
    'Test::TempDir::Tiny' => '0',
  },
  'configure_requires' => {
    'ExtUtils::CppGuess' => '0.07',
    'Module::Build' => '0.38',
    'Module::Build::WithXSpp' => '0.10'
  },
  'test_requires' => {
    'Test::TempDir::Tiny' => '0.016'
  },
  # HAS_BOOL         : stops Perl/lib/CORE/handy.h from doing "#  define bool char" for MSVC
  extra_compiler_flags => [qw(-DHAS_BOOL -I/usr/local/boost_1_60_0)],
  # Provides extra C typemaps that are auto-merged
  extra_typemap_modules => {
    'ExtUtils::Typemaps::Default' => '0.05',
  },
  #for MSVC builds
  early_includes => [qw(
    cstring
    cstdlib
    ostream
    )]
);

$build->create_build_script;

use Test::More;

require_ok 'Math::MatrixReal';
require_ok 'Algorithm::QuasiNewton';

my $f = sub {
   my $vector = shift;
   my $x = $vector->element(1,1);
   my $y = $vector->element(2,1);
   my $z = $vector->element(3,1);

   return ($x - 1) * ($x - 1) + ($y - 2) * ($y - 2) + ($z - 3) * ($z - 3);

};

my $df = sub {
   my $vector = shift;
   my $x = $vector->element(1,1);
   my $y = $vector->element(2,1);
   my $z = $vector->element(3,1);

   return Math::MatrixReal->new_from_cols([ [2 * $x - 2, 2 * $y - 4, 2 * $z - 6] ]);

};

my $init_param = Math::MatrixReal->new_from_cols([ [0,0,0] ]);
my $newton = Algorithm::QuasiNewton->new(f => $f, df => $df, x => $init_param, algorithm => 'lbfgs');
my $result = $newton->run();

is(sprintf(".6f",$result->element(1,1)),sprintf(".6f",1.0));
is(sprintf(".6f",$result->element(2,1)),sprintf(".6f",2.0));
is(sprintf(".6f",$result->element(3,1)),sprintf(".6f",3.0));

my $init_param = Math::MatrixReal->new_from_cols([ [0,0,0] ]);
my $newton = Algorithm::QuasiNewton->new(f => $f, df => $df, x => $init_param, algorithm => 'bfgs');
my $result = $newton->run();

is(sprintf(".6f",$result->element(1,1)),sprintf(".6f",1.0));
is(sprintf(".6f",$result->element(2,1)),sprintf(".6f",2.0));
is(sprintf(".6f",$result->element(3,1)),sprintf(".6f",3.0));

done_testing();


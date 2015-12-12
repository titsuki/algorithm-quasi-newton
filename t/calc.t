use Test::More;
use Data::Dumper;

require_ok 'Algorithm::QuasiNewton';

my $f = sub {
   my ($x,$y,$z) = @_;
   return $x * $x * $x + $y * $y + $z;
};

my $df = sub {
   my ($x,$y,$z) = @_;
   # x^3 + y^2 + z
   #d/dx -> 3x^2
   #d/dy -> 2y
   #d/dz -> 1
   return [3 * $x * $x, 2 * $y, 1];
};

my $newton = Algorithm::QuasiNewton->new(f => $f, df => $df, x => [0,0,0]);
my $result = $newton->run();
print STDERR Dumper($result);

done_testing();


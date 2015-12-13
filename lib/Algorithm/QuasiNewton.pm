package Algorithm::QuasiNewton;

use Mouse;
use Math::MatrixReal;

our $EPS = 1e-10;

has 'f' => (
    is => 'ro',
    isa => 'CodeRef'
    );

has 'df' => (
    is => 'ro',
    isa => 'CodeRef'
    );

has 'x' => (
    is => 'rw',
    isa => 'Math::MatrixReal'
    );

sub run {
    my $self = shift;
    my ($rows, $columns) = $self->{x}->dim();
    my $B = Math::MatrixReal->new_diag([map { 1;} @{ [1..$rows ] }]);

    my $fx = $self->f->($self->{x});
    my $g = $self->df->($self->{x});
    my $prev_fx;
    do {
	$prev_fx = $fx;
	my $gradient_direction = -1 * $B * $g;

	my $prev_x = $self->{x};
	$self->{x} = $self->golden_section_search($gradient_direction);
	
	$fx = $self->f->($self->{x});

	my $prev_g = $g->clone();
	$g = $self->df->($self->{x});

	my $y = $g - $prev_g;
	my $s = $self->{x} - $prev_x;

	my $ys = (~$y * $s)->element(1,1);
	my $yBy = (~$y * $B * $y)->element(1,1);
	$B = $B
	    + (1.0 + $yBy / $ys) * (1.0 / $ys) * ($s * ~$s)
	    - (1.0 / $ys) * ($B * $y * ~$s + $s * ~$y * ~$B);
    } while(abs($fx - $prev_fx) > $EPS);
    return $self->{x};
}

sub golden_section_search {
    my ($self,$gradient_direction) = @_;

    my $p = $self->{x};
    my $r = 2 / (3 + sqrt(5));
    my $a = 0;
    my $b = 1;
    my $t = $r * ($b - $a);
    my $c = $a + $t;
    my $d = $b - $t;
    my $fc;
    my $fd;
    
    $fc = $self->f->($p + $c * $gradient_direction);
    $fd = $self->f->($p + $d * $gradient_direction);
    
    while($d - $c > $EPS){
	if($fc > $fd){
	    my $a = $c;
	    $c = $d;
	    $d = $b - $r * ($b - $a);
	    $fc = $fd;
	    $fd = $self->f->($p + $d * $gradient_direction);
	} else {
	    my $b = $d;
	    $d = $c;
	    $c = $a - $r * ($b - $a);
	    $fd = $fc;
	    $fc = $self->f->($p + $c * $gradient_direction);
	}
    }
    return $p + $d * $gradient_direction;
}

1;

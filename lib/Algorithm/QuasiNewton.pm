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

has 'algorithm' => (
    is => 'rw',
    isa => 'Str',
    default => 'lbfgs'
    );

has 'max_iteration' => (
    is => 'ro',
    isa => 'Int',
    default => 1000
    );

has 'memory' => (
    is => 'ro',
    isa => 'Int',
    default => 10
    );

sub run {
    my $self = shift;
    if($self->{algorithm} eq 'bfgs'){
	return $self->bfgs();
    }
    else {
	return $self->lbfgs();
    }
}

sub lbfgs {
    my $self = shift;
    
    my ($rows, $columns) = $self->{x}->dim();
    my $B = Math::MatrixReal->new_diag([map { 1;} @{ [1..$rows ] }]);
    my $g = $self->df->($self->{x});
    my $fx = $self->f->($self->{x});

    for(my $iter = 0; $iter <= $self->{max_iteration}; $iter++){
	my $prev_fx = $fx;
	my $prev_x = $self->{x};
	my $prev_g = $g->clone();

	$g = $self->df->($self->{x});
	$fx = $self->f->($self->{x});

	if(abs($fx - $prev_fx) <= $EPS){
	    last;
	}

	my $bound = ($iter - $self->{memory} < 0 ? 0 : $iter - $self->{memory});
	my $q = $g;

	my $y = $g - $prev_g;
	my $s = $self->{x} - $prev_x;
	my $ys = ~$y * $s;

	my $rho = 0;
	if($ys->element(1,1) != 0){
	    $rho = 1.0 / $ys->element(1,1);
	}

	my $alpha;
	for(my $i = $iter - 1; $i >= $bound; $i--){
	    $alpha->{$i} = $rho * (~$s * $q)->element(1,1);
	    $q = $q - $y * $alpha->{$i};
	}

	my $r;

	$r = $q;
	for(my $i = $bound; $i < $iter; $i++){
	    my $beta = $rho * (~$y * $r)->element(1,1);	    
	    $r = $r + $s * ($alpha->{$i} - $beta);
	}

	my $gradient_direction = -$r;
	$self->{x} = $self->golden_section_search($gradient_direction);
    }
    return $self->{x};
}

sub bfgs {
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
	    $a = $c;
	    $c = $d;
	    $d = $b - $r * ($b - $a);
	    $fc = $fd;
	    $fd = $self->f->($p + $d * $gradient_direction);
	} else {
	    $b = $d;
	    $d = $c;
	    $c = $a + $r * ($b - $a);
	    $fd = $fc;
	    $fc = $self->f->($p + $c * $gradient_direction);
	}
    }
    return $p + $d * $gradient_direction;
}

1;

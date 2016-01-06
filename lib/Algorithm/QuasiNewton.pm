package Algorithm::QuasiNewton;

use Mouse;
use Math::MatrixReal;
use Data::Dumper;
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

has 'mode' => (
    is => 'rw',
    isa => 'Str',
    default => 'lbfgs'
    );

sub run {
    my $self = shift;
    if($self->{mode} eq 'bfgs'){
	return $self->bfgs();
    }
    else {
	return $self->lbfgs();
    }
}

sub lbfgs {
    my $self = shift;
    
    my $m = 100;
    my $incr;
    my $bound;
    my ($rows, $columns) = $self->{x}->dim();
    my $B = Math::MatrixReal->new_diag([map { 1;} @{ [1..$rows ] }]);
    my $fx = $self->f->($self->{x});
    my $g = $self->df->($self->{x});;
    my $prev_g = $g;
    my $prev_fx = $fx;
    my $prev_x = $self->{x};

    for(my $iter = 0; $iter <= 1000; $iter++){
	if($iter <= $m){
	    $incr = 0;
	    $bound = $iter;
	}
	else{
	    $incr = $iter - $m;
	    $bound = $m;
	}

	my $q = $g;
	my $alpha;
	for(my $i = $bound - 1; $i >= 0; $i--){
	    my $j = $i + $incr;
	    my $y = $g - $prev_g;
	    my $s = $self->{x} - $prev_x;
	    my $ys = ~$y * $s;
	    print STDERR $ys;
	    my $rho = ~$ys;
	    $alpha->{$bound} = ($rho * ~$s * $q)->element(1,1);
	    print STDERR Dumper($alpha->{$bound});
	    # print STDERR Dumper($q);
	    # print STDERR Dumper($y);
	    $q = $q - $alpha->{$bound} * $y;
	}

	my $r = $B * $q;
	for(my $i = 0; $i < $bound; $i++){
	    my $j = $i + $incr;
	    my $y = $g - $prev_g;
	    my $s = $self->{x} - $prev_x;
	    my $ys = ~$y * $s;
	    my $rho = ~$ys;
	    my $beta = ($rho * ~$y * $r)->element(1,1);
	    $r = $r + $s * ($alpha->{$i} - $beta);
	}
	my $gradient_direction = $r;
	$prev_fx = $fx;
	$self->{x} = $self->golden_section_search($gradient_direction);
	
	$fx = $self->f->($self->{x});

	$prev_g = $g->clone();
	$g = $self->df->($self->{x});
    }
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

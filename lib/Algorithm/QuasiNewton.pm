package Algorithm::QuasiNewton;

use Mouse;

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
    isa => 'ArrayRef[Num]'
    );

sub run {
    my $self = shift;
    my $B;
    
    for(my $i = 0; $i < @{ $self->{x} }; $i++){
	for(my $j = 0; $j < @{ $self->{x} }; $j++){
	    $B->[$i]->[$j] = 0;
	}
    }

    for(my $i = 0; $i < @{ $self->{x} }; $i++){
	$B->[$i]->[$i] = 1;
    }

    my $fx = $self->f->(@{ $self->{x} });
    my $g = $self->df->(@{ $self->{x} });
    my $s;
    my $fz;
    do {
	$fz = $fx;
	for(my $i = 0; $i < @{ $self->{x} }; $i++){
	    for(my $j = 0; $j < @{ $self->{x} }; $j++){
		$s->[$i] -= $B->[$i]->[$j] * $g->[$j];
	    }
	}

	$self->golden_section_search($s);
	
	$fx = $self->f->($self->{x});
	my $y = $g;
	$g = $self->df->($self->{x});
	for(my $i = 0; $i < @{ $self->{x} }; $i++){
	    $y->[$i] = $g->[$i] - $y->[$i];
	}

	my $By;
	for(my $i = 0; $i < @{ $self->{x} }; $i++){
	    for(my $j = 0; $j < @{ $self->{x} }; $j++){
		$By->[$i] += $B->[$i]->[$j] * $y->[$j];
	    }
	}
	my $sy = $self->dot($s,$y) + $EPS;
	my $yBy = $self->dot($y,$By) + $EPS;
	my $u;
	for(my $i = 0; $i < @{ $By }; $i++){
	    $u->[$i] = $s->[$i]/$sy - $By->[$i]/$yBy;
	}

	for(my $i = 0; $i < @{ $s }; $i++){
	    for(my $j = 0; $j < @{ $s }; $j++){
		$B->[$i]->[$j] += $s->[$i] * $s->[$j]/$sy
		    - $By->[$i] * $By->[$j]/$yBy
		    + $u->[$i] * $u->[$j] * $yBy;
	    }
	}
    } while(abs($fx - $fz) > $EPS);
    return $self->{x};
}

sub golden_section_search {
    my ($self,$s) = @_;

    my $p = $self->{x};
    my $r = 2 / (3 + sqrt(5));
    my $a = 0;
    my $b = 1;
    my $t = $r * ($b - $a);
    my $c = $a + $t;
    my $d = $b - $t;
    my $fc;
    my $fd;
    # a c d b
    for(my $i = 0; $i < @{ $self->{x} }; $i++){
	$self->{x}->[$i] = $p->[$i] + $c * $s->[$i];
    }
    $fc = $self->f->(@{ $self->{x} });
    
    for(my $i = 0; $i < @{ $self->{x} }; $i++){
	$self->{x}->[$i] = $p->[$i] + $d * $s->[$i];
    }
    $fd = $self->f->(@{ $self->{x} });
    
    while($d - $c > $EPS){
	if($fc > $fd){
	    my $a = $c;
	    $c = $d;
	    $d = $b - $r * ($b - $a);
	    for(my $i = 0; $i < @{ $self->{x} }; $i++){
		$self->{x}->[$i] = $p->[$i] + $d * $s->[$i];
	    }
	    $fc = $fd;
	    $fd = $self->f->($self->{x});
	} else {
	    my $b = $d;
	    $d = $c;
	    $c = $a - $r * ($b - $a);
	    for(my $i = 0; $i < @{ $self->{x} }; $i++){
		$self->{x}->[$i] = $p->[$i] + $c * $s->[$i];
	    }
	    $fd = $fc;
	    $fc = $self->f->($self->{x});
	}
    }
}

sub dot {
    my($self,$x,$y) = @_;
    my $res = 0;
    for(my $i = 0; $i < @{ $self->{x} }; $i++){
	$res += $x->[$i] * $y->[$i];
    }
    return $res;
}

1;

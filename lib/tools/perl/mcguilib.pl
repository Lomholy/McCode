use Tk;
use Tk::DialogBox;
use Tk::ROText;
use Tk::Listbox;

sub get_dir_name {
    my ($dlg, $default) = @_;
    my $oldgrab = $dlg->grabStatus;
    $dlg->grabRelease;
    my $f = $default ?
	$dlg->getSaveFile(-title => "Select output file name",
			  -initialfile => $default) :
	$dlg->getSaveFile(-title => "Select output file name");
    $dlg->grab if $oldgrab eq 'local';
    return $f;
}

# Query user for instrument parameters and simulation options for a
# McStas simulation.
# Input: top-level window for the dialog, instrument info descriptor, and
#        simulation info descriptor.
# Output: user action ("Start" or "Cancel") and new simulation info
#         descriptor.

my %typeabbrev = ('double' => "D", 'int' => "I", 'string' => "S");
my $typehelp = "(D=floating point, I=integer, S=string)";

sub simulation_dialog {
    my ($win, $ii, $origsi) = @_;
    my %si = $origsi ? %$origsi : ();
    my $doseed;
    if($origsi->{'Seed'}) {
	$si{'Seed'} = $origsi->{'Seed'};
	$doseed = 1;
    } else {
	$si{'Seed'} = "";
	$doseed = 0;
    }
    $si{'Autoplot'} = 0 unless $si{'Autoplot'};
    $si{'Ncount'} = 1e6 unless $si{'Ncount'};
    $si{'Trace'} = 0 unless $si{'Trace'};

    my $dlg = $win->DialogBox(-title => "Run simulation",
			      -buttons => ["Start", "Cancel"]);

    $dlg->add('Label',
	      -text => "Instrument source: $ii->{'Instrument-source'}",
	      -anchor => 'w',
	      -justify => 'left')->pack(-fill => 'x');

    # Set up the parameter input fields.
    my @parms = @{$ii->{'Parameters'}};
    my $numrows = int ((@parms + 2)/3);
    if($numrows > 0) {
	$dlg->add('Label',
		  -text => "Instrument parameters $typehelp:",
		  -anchor => 'w',
		  -justify => 'left')->pack(-fill => 'x');
	my $parm_frame = $dlg->add('Frame');
	$parm_frame->pack;
	my $row = 0;
	my $col = 0;
	my $p;
	for $p (@parms) {
	    # Give parameter type as abbrevation.
	    my ($type, $u);
	    $type = $ii->{'Parameter-types'}{$p};
	    $u = " ($typeabbrev{$type})" if $type;
	    $u = "" unless $u;
	    my $w = $parm_frame->Label(-text => "$p$u:", -justify => 'right');
	    $w->grid(-row => $row, -column => $col, -sticky => 'e');
	    $col++;
	    $si{'Params'}{$p} = "" unless defined($si{'Params'}{$p});
	    $w = $parm_frame->Entry(-relief => 'sunken',
				    -width=>10,
				    -textvariable => \$si{'Params'}{$p},
				    -justify => 'right');
	    $w->grid(-row => $row, -column => $col);
	    $col++;
	    if($col >= 6) {
		$col = 0;
		$row++;
	    }
	}
    }

    my $opt_frame = $dlg->add('Frame');
    $opt_frame->pack(-anchor => 'w');

    my $f0 = $opt_frame->Frame;
    $f0->pack(-anchor => 'w');
    $f0->Label(-text => "Output to:")->pack(-side => 'left');
    my $dir_entry = $f0->Entry(-relief => 'sunken',
			       -width=>30,
			       -justify => 'left',
			       -textvariable => \$si{'Dir'});
    $dir_entry->pack(-side => 'left');
    $f0->Button(-text => "Browse ...",
		-command => sub { my $d = get_dir_name($dlg, $si{'Dir'});
				  $si{'Dir'} = $d if $d; } )->pack;

    my $f1 = $opt_frame->Frame;
    $f1->pack(-anchor => 'w');
    $f1->Label(-text => "Neutron count:")->pack(-side => 'left');
    my $ncount_entry = $f1->Entry(-relief => 'sunken',
				  -width=>10,
				  -justify => 'right',
				  -textvariable => \$si{'Ncount'});
    $ncount_entry->pack(-side => 'left');
    $opt_frame->Checkbutton(-text => "Plot results",
			    -variable => \$si{'Autoplot'},
			    -relief => 'flat')->pack(-anchor => 'w');
    my $f2 = $opt_frame->Frame;
    $f2->pack(-anchor => 'w');
    $f2->Radiobutton(-text => "Random seed",
		     -variable => \$doseed,
		     -relief => 'flat',
		     -value => 0)->pack(-side => 'left');
    $f2->Radiobutton(-text => "Set seed to",
		     -variable => \$doseed,
		     -relief => 'flat',
		     -value => 1)->pack(-side => 'left');
    $ncount_entry = $f2->Entry(-relief => 'sunken',
			       -width=>10,
			       -textvariable => \$si{'Seed'},
			       -justify => 'right');
    $ncount_entry->pack(-side => 'left');

    my $f3 = $opt_frame->Frame;
    $f3->pack(-anchor => 'w');
    $f3->Radiobutton(-text => "Simulate",
		     -variable => \$si{'Trace'},
		     -relief => 'flat',
		     -value => 0)->pack(-side => 'left');
    $f3->Radiobutton(-text => "Trace",
		     -variable => \$si{'Trace'},
		     -relief => 'flat',
		     -value => 1)->pack(-side => 'left');

    my $res = $dlg->Show;
    $si{'Seed'} = 0 unless $doseed;
    return ($res, \%si);
}

my $current_plot;

sub dialog_plot_single {
    my ($cl,$di) = @_;
    $current_plot = $cl->index('active');
    single_plot("/xserv", $di->[$current_plot], 0);
}

sub dialog_hardcopy {
    my ($dlg, $di, $type) = @_;
    my $default = $current_plot == -1 ?
	"mcstas.ps" :
	($di->[$current_plot]{'Component'} . ".ps");
    my $oldgrab = $dlg->grabStatus;
    $dlg->grabRelease;
    my $f = $dlg->getSaveFile(-defaultextension => "ps",
			      -title => "Select postscript file name",
			      -initialfile => $default);
    $dlg->grab if $oldgrab eq 'local';
    return 0 unless $f;
    if($current_plot == -1) {
	overview_plot("\"$f\"/$type", $di, 0);
    } else {
	my $comp = $di->[$current_plot]{'Component'};
	single_plot("\"$f\"/$type", $di->[$current_plot], 0);
    }
}
    
sub plot_dialog {
    my ($win, $ii, $si, $di, $sim_file_name) = @_;
    $current_plot = -1;	# Component index, or -1 -> overview.
    my $dlg = $win->DialogBox(-title => "Plot results",
			      -buttons => ["Close"]);

    my $lf = $dlg->add('Frame');
    $lf->pack(-side => 'left');
    $lf->Label(-text => "Monitors and detectors",
	       -anchor => 'w')->pack(-fill => 'x');
    my $cl = $lf->Scrolled('Listbox',
			  -width => 25,
			  -height => 10,
			  -setgrid => 1,
			  -scrollbars => 'se');
    $cl->pack(-expand => 'yes', -fill => 'y', -anchor => 'w');
    $cl->bind('<Double-Button-1>' => sub { dialog_plot_single($cl,$di);
				       $dlg->raise; } );
    $cl->insert(0, map $_->{'Component'}, @$di);
    $cl->activate(0);
    $lf->Button(-text => "Plot",
		-command => sub { dialog_plot_single($cl,$di);
			          $dlg->raise; } )->pack;
    $lf->Button(-text => "Overview plot",
		-command => sub {
		    overview_plot("/xserv", $di, 0);
		    $dlg->raise;
		    $current_plot = -1; }
		)->pack;
    $lf->Button(-text => "B&W postscript",
		-command => sub { dialog_hardcopy($dlg,
						  $di, "ps"); }
		)->pack;
    $lf->Button(-text => "Colour postscript",
		-command => sub { dialog_hardcopy($dlg,
						  $di, "cps"); }
		)->pack;
#     $lf->Button(-text => "Select from overview",
# 		-command => sub {
# 		    my ($c, $idx) = overview_plot("/xserv", $di, 1);
# 		    $cl->activate($idx);
# 		    $current_plot = -1;}
# 		)->pack;
    my $rf = $dlg->add('Frame');
    $rf->pack(-side => 'top');
    $rf->Label(-text => <<END,
Date: $si->{'Date'}
Instrument name: $ii->{'Name'}
Source: $ii->{'Instrument-source'}
Neutron count: $si->{'Ncount'}
Simulation file: $sim_file_name
END
	       -anchor => 'w',
	       -justify => 'left')->pack(-fill => 'x');

    overview_plot("/xserv", $di, 0);
    my $res = $dlg->Show;
    return ($res);
}

sub fetch_comp_info {
    my ($cname, $cinfo) = @_;
    unless($cinfo->{$cname}) {
	my $info = component_information($cname);
	return undef unless $info;
	$cinfo->{$cname} = $info;
    }
    return $cinfo->{$cname};
}

sub comp_instance_dialog {
    my ($w, $comp) = @_;
    my ($i, $j, $p);
    my $r = {
	'INSTANCE' => "",
	'DEFINITION' => "",
	'VALUE' => { },
	'AT' => { 'x' => "", 'y' => "", 'z' => "", 'relative' => "" },
	'ROTATED' => { 'x' => "", 'y' => "", 'z' => "", 'relative' => "" }
    };
    my $dlg = $w->Toplevel(-title => "$comp->{'name'}");
#    $dlg->transient($dlg->Parent->toplevel);
    $dlg->withdraw;
    # Add labels
    $dlg->Label(-text => "Component definition: $comp->{'name'}",
		-anchor => 'w')->pack(-fill => 'x');
    $r->{'DEFINITION'} = $comp->{'name'};
    $dlg->Label(-text => "$comp->{'identification'}{'short'}",
		-padx => 32, -anchor => 'w')->pack(-fill => 'x');
    my $f1 = $dlg->Frame();
    $f1->pack(-fill => 'x');
    $f1->Label(-text => "Author: $comp->{'identification'}{'author'}",
	       -anchor => 'w', -justify => 'left')->pack(-side => 'left');
    $f1->Label(-text => "Date: $comp->{'identification'}{'date'}",
	       -anchor => 'w', -justify => 'left')->pack(-side => 'right');
    $dlg->Label(-text => "Origin: $comp->{'identification'}{'origin'}",
		-anchor => 'w', -justify => 'left')->pack(-fill => 'x');
    my $f2 = $dlg->Frame();
    $f2->pack(-fill => 'x');
    $f2->Label(-text => "Instance name: ", -fg => 'blue',
	       -anchor => 'w', -justify => 'left')->pack(-side => 'left');
    my $entry = $f2->Entry(-relief => 'sunken', -width => 20,
			   -textvariable => \$r->{'INSTANCE'},
			   -justify => 'left');
    $entry->pack(-side => 'left');
    $entry->focus;
    my $t = $dlg->Scrolled(qw/ROText -relief sunken -bd 2 -setgrid true
			   -height 18 -width 80 -scrollbars osoe/);
    $t->pack(-expand => 'yes', -fill => 'both');
    $t->tagConfigure('BLUE', -foreground => 'blue');
    $t->tagConfigure('RED', -foreground => 'red');
    $t->insert('end', "PARAMETERS:\n\n", 'RED');
    for $p (@{$comp->{'inputpar'}}) {
	$t->insert('end', "$p:", 'BLUE');
	my $entry = $t->Entry(-relief => 'sunken', -width => 10,
			      -textvariable => \$r->{'VALUE'}{$p},
			      -justify => 'right');
	$t->window('create', 'end', -window => $entry);
	my $unit = $comp->{'parhelp'}{$p}{'unit'};
	if($unit) {
	    $t->insert('end', " [");
	    $t->insert('end', $unit, 'RED');
	    $t->insert('end', "]");
	}
	my $def = $comp->{'parhelp'}{$p}{'default'};
	$t->insert('end', " (OPTIONAL, default $def)") if defined($def);
	$t->insert('end', "\n");
	$t->insert('end', $comp->{'parhelp'}{$p}{'text'})
	    if $comp->{'parhelp'}{$p}{'text'};
	$t->insert('end', "\n\n");
    }
    $t->insert('end', "DESCRIPTION:\n\n", 'RED');
    $t->insert('end', $comp->{'description'});
    $t->see("1.0");
    $t->markSet('insert', "1.0");

    my %delim = ('x' => ", ", 'y' => ", ", 'z' => "");
    for $j (qw/AT ROTATED/) {
	my $f3 = $dlg->Frame();
	$f3->pack(-fill => 'x');
	$f3->Label(-text => "$j", -fg => 'blue',
		   -anchor => 'w', -justify => 'left')->pack(-side => 'left');
	$f3->Label(-text => " (",
		   -anchor => 'w', -justify => 'left')->pack(-side => 'left');
	for $i (qw/x y z/) {
	    my $entry = $f3->Entry(-relief => 'sunken', -width => 6,
				   -textvariable => \$r->{$j}{$i},
				   -justify => 'right');
	    $entry->pack(-side => 'left');
	    $f3->Label(-text => $delim{$i})->pack(-side => 'left');
	}
	$f3->Label(-text => ")  ")->pack(-side => 'left');
	$f3->Label(-text => "RELATIVE ",
		   -fg => 'blue')->pack(-side => 'left');
	my $entry2 = $f3->Entry(-relief => 'sunken', -width => 12,
				-textvariable => \$r->{$j}{'relative'},
				-justify => 'right');
	$entry2->pack(-side => 'left');
    }

    my $bot_frame = $dlg->Frame(-relief => "raised", -bd => 1);
    $bot_frame->pack(-side => "top", -fill => "both",
		     -ipady => 3, -ipadx => 3);
    my $selected;
    my $ok_cmd = sub { $selected = 'OK' };
    my $cancel_cmd = sub { $selected = 'CANCEL' };
    my $okbut = $bot_frame->Button(-text => "Ok", -command => $ok_cmd);
    $okbut->pack(-side => "left", -expand => 1, -padx => 1, -pady => 1);
    my $cancelbut = $bot_frame->Button(-text => "Cancel",
				       -command => $cancel_cmd);
    $cancelbut->pack(-side => "left", -expand => 1, -padx => 1, -pady => 1);
    $dlg->protocol("WM_DELETE_WINDOW" => $cancel_cmd);
    $dlg->bind('<Escape>' => $cancel_cmd);
    $dlg->bind('<Return>' => $ok_cmd);

    my $old_focus = $dlg->focusSave;
    my $old_grab = $dlg->grabSave;
    $dlg->Popup;
    $dlg->grab;
    $dlg->waitVariable(\$selected);
    $dlg->grabRelease;
    $dlg->destroy;
    &$old_focus;
    &$old_grab;
    return ($selected eq 'OK' ? $r : undef);
}

# Make a special Listbox class that invokes a user-specified callback
# each time an item is selected.
{
    package Tk::MyListbox;
    @Tk::MyListbox::ISA = qw/Tk::Listbox/; # Inherit from Tk::Listbox.
    Tk::Widget->Construct('MyListbox');
    sub ClassInit {
	my($cw, @args) = @_;
	$cw->SUPER::ClassInit(@args);
    }
    sub Populate {
	my($cw, @args) = @_;
	$cw->SUPER::Populate(@args);
    }
    my $selecthook;
    # Set the callback to invoke whenever an item is selected.
    sub selecthook {
	my ($w, $f) = @_;
	$selecthook = $f;
    }
    # Override the selectionSet() method, which gets called every time
    # a Listbox selection is made. This seems to be a private method,
    # hopefully it will not change in a later PerlTk version.
    sub selectionSet {
	my ($w,@args) = @_;
	my @r = $w->SUPER::selectionSet(@args);
	&$selecthook() if $selecthook; # Invoke user callback
	return @r;
    }
}

sub comp_select_dialog {
    my ($w, $clist, $cinfo) = @_;
    my $dlg = $w->Toplevel(-title => "Select component definition");
    $dlg->transient($dlg->Parent->toplevel);
    $dlg->withdraw;
    my $f = $dlg->Frame();
    $f->pack(-side => 'top');
    $f->Label(-text => "Available component definitions:")->pack;
    my $list = $f->Scrolled('MyListbox', -width => 50, -height => 10,
			    -setgrid => 1, -scrollbars => 'osre');
    $list->pack(-expand => 'yes', -fill => 'y', -anchor => 'n');
    my @sorted = sort {compname($a) cmp compname($b)} @$clist;
    my @namelist = map compname($_), @sorted;
    $list->insert(0, @namelist);
    $list->activate(0);
    my $name = $f->Label(-text => "Name: ", -anchor => 'w');
    $name->pack(-fill => 'x');
    my $loc = $f->Label(-text => "Location: ", -anchor => 'w');
    $loc->pack(-fill => 'x');
    my $text = $f->Scrolled(qw/ROText -relief sunken -bd 2 -setgrid true
			    -height 10 -width 80 -scrollbars osoe/);
    $text->pack();
    $text->tagConfigure('SHORT', -foreground => 'blue');
    my $f1 = $f->Frame();
    $f1->pack(-fill => 'x');
    my $author = $f1->Label(-text => "Author: ",
			    -anchor => 'w', -justify => 'left');
    $author->pack(-side => 'left');
    my $date = $f1->Label(-text => "Date: ",
			  -anchor => 'w', -justify => 'left');
    $date->pack(-side => 'right');

    my $bot_frame = $dlg->Frame(-relief => "raised", -bd => 1);
    $bot_frame->pack(-side => "top", -fill => "both",
		     -ipady => 3, -ipadx => 3);
    my $selected;
    my $select_cmd = sub {
	my $cname = $sorted[$list->curselection()];
	my $info = fetch_comp_info($cname, $cinfo);
	$name->configure(-text => "Name: $info->{'name'}");
	$loc->configure(-text => "Location: $cname");
	$author->configure(-text =>
			   "Author: $info->{'identification'}{'author'}");
	$date->configure(-text => "Date: $info->{'identification'}{'date'}");
	$text->delete("1.0", "end");
	$text->insert("end", "$info->{'identification'}{'short'}\n\n","short");
	$text->insert("end", $info->{'description'});
    };
    my $accept_cmd = sub { $selected = 'Ok'; };
    my $cancel_cmd = sub { $selected = 'Cancel'; };
    # Set up function to be called each time an item is selected.
    $list->selecthook($select_cmd);
    $list->bind('<Double-Button-1>' => $accept_cmd);
    $list->bind('<Return>' => $accept_cmd);
    my $okbut = $bot_frame->Button(-text => "Ok", -command => $accept_cmd);
    $okbut->pack(-side => "left", -expand => 1, -padx => 1, -pady => 1);
    my $cancelbut = $bot_frame->Button(-text => "Cancel",
				       -command => $cancel_cmd);
    $cancelbut->pack(-side => "left", -expand => 1, -padx => 1, -pady => 1);
    $dlg->protocol("WM_DELETE_WINDOW" => $cancel_cmd);
    $dlg->bind('<Escape>' => $cancel_cmd);
    $dlg->bind('<Return>' => $accept_cmd);

    my $old_focus = $dlg->focusSave;
    my $old_grab = $dlg->grabSave;
    $list->focus;
    $dlg->Popup;
    $dlg->grab;
    $dlg->waitVariable(\$selected);
    my $selected_comp = ($selected eq 'Ok' ?
			 $sorted[$list->curselection()] :
			 undef);
    $dlg->grabRelease;
    $dlg->destroy;
    &$old_focus;
    &$old_grab;
    return $selected_comp;
}

1;

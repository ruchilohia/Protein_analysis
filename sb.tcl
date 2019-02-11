proc sb { molid dir } {
     
     source /u1/home/lohia/data/utilities/vis/sscache.tcl
     start_sscache $molid
     set sel_salt [atomselect $molid "protein"]
     atomselect macro acidic "resname ASP GLU"
     atomselect macro basic "resname ARG LYS"
     package require saltbr
     saltbr -sel $sel_salt -upsel yes -ondist 100 -writefiles yes -outdir $dir -log sb_log -frames 1:end
     

 } 
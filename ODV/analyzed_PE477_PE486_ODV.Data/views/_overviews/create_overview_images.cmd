# Creates full-range overview plots (collection assumed open).

set_var COLLDIR  = C:/Users/hisham.shaikh/OneDrive - UGent/Projects/PE477_PE486_Viral_Production/ODV
set_var COLLNAME = analyzed_PE477_PE486_ODV
set_var VIEWSDIR = %COLLDIR%/%COLLNAME%.Data/views/_overviews

load_view %VIEWSDIR%/OVERVIEW_AllStationsMap.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_AllStationsMap.gif, 300

load_view %VIEWSDIR%/OVERVIEW_vars_0002-0028.xview
set_axis_full_ranges -1
set_axis_full_ranges 0
save_view %VIEWSDIR%/OVERVIEW_vars_0002-0028.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_vars_0002-0028.gif, 300

load_view %VIEWSDIR%/OVERVIEW_vars_0029-0055.xview
set_axis_full_ranges -1
set_axis_full_ranges 0
save_view %VIEWSDIR%/OVERVIEW_vars_0029-0055.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_vars_0029-0055.gif, 300

load_view %VIEWSDIR%/OVERVIEW_vars_0056-0082.xview
set_axis_full_ranges -1
set_axis_full_ranges 0
save_view %VIEWSDIR%/OVERVIEW_vars_0056-0082.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_vars_0056-0082.gif, 300

load_view %VIEWSDIR%/OVERVIEW_vars_0083-0109.xview
set_axis_full_ranges -1
set_axis_full_ranges 0
save_view %VIEWSDIR%/OVERVIEW_vars_0083-0109.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_vars_0083-0109.gif, 300

load_view %VIEWSDIR%/OVERVIEW_vars_0110-0136.xview
set_axis_full_ranges -1
set_axis_full_ranges 0
save_view %VIEWSDIR%/OVERVIEW_vars_0110-0136.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_vars_0110-0136.gif, 300

load_view %VIEWSDIR%/OVERVIEW_vars_0137-0163.xview
set_axis_full_ranges -1
set_axis_full_ranges 0
save_view %VIEWSDIR%/OVERVIEW_vars_0137-0163.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_vars_0137-0163.gif, 300

load_view %VIEWSDIR%/OVERVIEW_vars_0164-0173.xview
set_axis_full_ranges -1
set_axis_full_ranges 0
save_view %VIEWSDIR%/OVERVIEW_vars_0164-0173.xview
export_graphics -1, %COLLDIR%/%COLLNAME%_OVERVIEW_vars_0164-0173.gif, 300


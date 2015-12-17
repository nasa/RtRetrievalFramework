#!/usr/bin/env python

from types import ListType
from optparse import OptionParser
import L2_Input
import shutil
import copy
import re
import os

source_inp_files = ["oco_l2.inp", "oco_l2.run", "oco_l2.win"]
dst_inp_name = "oco_l2.inp"
dst_run_name = "oco_l2.run"
dst_win_name = "oco_l2.win"
bak_ext = ".bak_convert"

new_window_id = ['ABO2_Full', 'WCO2_Full', 'SCO2_Full']

def Find_Section(srch_objs, section_name):
    found_sects = []
    for inpt_obj in srch_objs:
        (names, sects) = inpt_obj.rootNode.Get_All_Sections()

        for (node_name, node_obj) in zip(names, sects):
            if node_name.lower() == section_name.lower():
                found_sects.append(node_obj)
    return found_sects

def Find_Keyword(srch_objs, keyword_name, default=None, error_if_unfound=True):
    found_keywords = []
    for inpt_obj in srch_objs:
        #        for key_obj in inpt_obj.rootNode.Get_All_Keyword_Objs():
        for sect_obj in inpt_obj.children:
            if sect_obj.type == 'assignment' and sect_obj.leaf.lower() == keyword_name.lower():
                found_keywords.append(sect_obj)
            elif sect_obj.type == 'section':
                for fnd_keys in Find_Keyword([sect_obj], keyword_name, None, False):
                    found_keywords.append(fnd_keys)
                    
    if default == None and error_if_unfound and len(found_keywords) == 0:
        raise ValueError('Could not find keyword: %s' % keyword_name)
    elif default != None and len(found_keywords) == 0:
        val_node = L2_Input.Node('value', leaf=default)
        key_node = L2_Input.Node('assignment', leaf=keyword_name, children=[val_node])
        found_keywords.append(key_node)
        
    return found_keywords

def Clean_Extra_Space(sect_obj):
    if not type(sect_obj) is ListType:
        sect_obj = [sect_obj]

    for curr_sect in sect_obj:
        for es_idx in range(len(curr_sect.endspace)):
            if curr_sect.endspace[es_idx].count('\n') > 1:
                curr_sect.endspace[es_idx] = curr_sect.endspace[es_idx].replace('\n', '') + '\n'
        
        for child in curr_sect.children:
            if child.type == 'assignment':
                subchilds = child.children
                if len(subchilds) > 0:
                    last_sub = subchilds[len(subchilds)-1]
                    if last_sub.endspace.count('\n') > 1:
                        last_sub.endspace = last_sub.endspace.replace('\n', '') + '\n'

def Clean_Assignment_Spacing(sect_obj):
    if not type(sect_obj) is ListType:
        sect_obj = [sect_obj]

    for curr_sect in sect_obj:
        max_key_len = 0
        for child in curr_sect.children:
            if child.type == 'assignment':
                max_key_len = max(max_key_len, len(child.leaf))
        for child in curr_sect.children:
            if child.type == 'assignment':
                child.endspace = ' ' * (max_key_len - len(child.leaf) + 1)
            
def Backup_Write(fileobj, filename):
    bak_filename = filename + bak_ext
    if os.path.exists(bak_filename):
        print 'Will not overwrite backup file: "%s"' % bak_filename
        print 'Skipping file with existing changes: "%s"' % filename
        return

    if os.path.exists(filename):
        shutil.copyfile(filename, bak_filename)

    fileobj.Write(filename, doIndent=True)


parser = OptionParser(usage="usage: %prog [options] <run_dir> [ <run_dir> ... ]")

parser.add_option( "-r", "--run_dir_file", dest="run_dir_file",
                   metavar="FILE",
                   help="file to read list of run directories from")

parser.add_option( "-w", "--update_win", dest="update_win",
                   default=False,
                   action="store_true",
                   help="Update %s windows file" % dst_win_name
                   )


# Parse command line arguments
(options, args) = parser.parse_args()


# Gather list of run directories to gather information from
run_dirs = []

if (len(args) > 0):
    for arg_dir in args:
        run_dirs.append(arg_dir)

if options.run_dir_file != None:
    if not os.path.exists(options.run_dir_file):
        parser.error("Run directory file '%s' does not exist" % options.run_dir_file)

    run_dir_fh = open(options.run_dir_file, 'r')
    for file_dir in run_dir_fh.readlines():
        run_dirs.append(file_dir.strip())
    run_dir_fh.close()

# Sort items from run dir list
run_dirs.sort()

if len(run_dirs) == 0:
    parser.error('at least one run directory must be specified for conversion')

for mod_dir in run_dirs:
    # Load source objects for searching
    source_objs = []
    for src_basename in source_inp_files:
        src_filename = "%s/%s" % (mod_dir, src_basename)
        if not os.path.exists(src_filename):
            raise IOError('Required file "%s" does not exist at "%s"' % (src_basename, mod_dir))
        source_objs.append( L2_Input.Input_File(src_filename) )

    ### generate new inp structure:
        
    dst_inp_fileobj = L2_Input.Input_File()

    # Create new window_info section
    win_info_sect = L2_Input.Section("WINDOW_INFO")
    win_info_sect.Add_Child(Find_Keyword(source_objs, 'spectral_window_file', 'oco_l2.win'))
    win_info_sect.Add_Child(Find_Keyword(source_objs, 'target_species', 'CO2'))
    spectral_windows_objs = Find_Keyword(source_objs, 'spectral_windows', error_if_unfound=False)
    win_count = 1
    for sw_obj in spectral_windows_objs:
        sw_obj.leaf = "%s(%d)" % (sw_obj.leaf, win_count)

        win_id = sw_obj.children[0].leaf
        if win_id.isdigit() and int(win_id) < len(new_window_id)+1:
            sw_obj.children = [ L2_Input.Node('value', new_window_id[int(win_id)-1]) ]
        
        win_info_sect.Add_Child(sw_obj)
        win_count += 1

    Clean_Extra_Space(win_info_sect)
    Clean_Assignment_Spacing(win_info_sect)   
    dst_inp_fileobj.rootNode.Add_Child(win_info_sect)

    # Find and copy existing sounding info section and add new elements
    sound_info_sect = Find_Section(source_objs, 'SOUNDING_INFO')[0]
    sound_info_sect.Add_Child(Find_Keyword(source_objs, 'noise_file', error_if_unfound=False))
    sound_info_sect.Add_Child(Find_Keyword(source_objs, 'pressure_file'))
    sound_info_sect.Add_Child(Find_Keyword(source_objs, 'absco_path'))
    sound_info_sect.Add_Child(Find_Keyword(source_objs, 'instrument', 'OCO'))
    Clean_Extra_Space(sound_info_sect)
    Clean_Assignment_Spacing(sound_info_sect)   
    dst_inp_fileobj.rootNode.Add_Child(sound_info_sect)

    # Find existing parameter definition section, fix retrievial indicies
    param_def_sect = Find_Section(source_objs, 'PARAMETER_DEFINITION')[0]
    param_def_sub_sects = param_def_sect.Get_All_Section_Nodes()
    for pd_sub_sect in param_def_sub_sects:
        keyword_objs = pd_sub_sect.Get_All_Keyword_Objs()
        ret_ind_objs = []
        for key_obj in keyword_objs:
            if key_obj.leaf == 'retrieval_indices':
                ret_ind_objs.append(key_obj)

        ret_ind_count = 1
        if len(ret_ind_objs) > 1:
            for ri_obj in ret_ind_objs:
                ri_obj.leaf =  "%s(%d)" % (ri_obj.leaf, ret_ind_count)
                ret_ind_count += 1

        Clean_Assignment_Spacing(pd_sub_sect)

    # This bit of code does the following:
    #   * finds ils files
    #   * makes a copy of function_type nodes 
    #     and delete originals from param block
    dst_ils_ap_name=''
    dst_ils_cov_name=''
    dst_ils_pert_name=''
    for pd_sub_sect in param_def_sub_sects:
        if pd_sub_sect.leaf[0].lower() == 'instrument':
            keyword_objs = pd_sub_sect.Get_All_Keyword_Objs()
            is_ils_type = False
            is_function_type = False
            for key_obj in keyword_objs:
                if key_obj.leaf.lower() == 'function_type':
                    nd_function_type = key_obj
                    pd_sub_sect.Del_Child(key_obj)
                elif key_obj.leaf.lower() == 'a_priori':
                    ils_ap_candidate = key_obj.children[0].leaf
                elif key_obj.leaf.lower() == 'covariance':
                    ils_cov_candidate = key_obj.children[0].leaf
                elif key_obj.leaf.lower() == 'perturb':
                    ils_pert_candidate = key_obj.children[0].leaf
                elif key_obj.leaf.lower() == 'type' \
                        and key_obj.children[0].leaf.lower() == 'ils':
                    is_ils_type = True
            if is_ils_type == True:
                dst_ils_ap_name = ils_ap_candidate
                dst_ils_cov_name = ils_cov_candidate
                dst_ils_pert_name = ils_pert_candidate
            else:
                dst_ils_ap_name =''
                dst_ils_cov_name = ''
                dst_ils_pert_name = ''

    if dst_ils_ap_name == '' and dst_ils_cov_name == '' and dst_ils_pert_name == '':
        print "No ILS block found; skip ILS conversion"
    elif dst_ils_ap_name == '' or dst_ils_cov_name == '' or dst_ils_pert_name == '':
        print "ERROR:  Not all ILS files could be resolved!"
        print "apriori=" + dst_ils_ap_name
        print "covariance=" + dst_ils_cov_name
        print "perturbation=" + dst_ils_pert_name
        os._exit(1)
    else:

        for dst_ils_curr_name in [ dst_ils_ap_name, dst_ils_cov_name, dst_ils_pert_name]:
            #########################################################
            if os.path.exists(dst_ils_curr_name):
                print "Updating '" + dst_ils_curr_name + "'"
            else:
                print "Skipping missing: %s" % dst_ils_curr_name
                continue
            #########################################################

            # write header block to 'ils.dat'
            dst_ils_fileobj = L2_Input.Input_File(dst_ils_curr_name)

            ils_source_objs = copy.copy(source_objs)
            ils_source_objs.append(dst_ils_fileobj)
            
            header_sect_srch = Find_Section([dst_ils_fileobj], 'HEADER')
            if header_sect_srch == None or len(header_sect_srch) == 0:
                raise IOError('Could not find section HEADER in file: %s' % dst_ils_curr_name)
            header_sect = header_sect_srch[0]
            header_sect.Add_Child(Find_Keyword(ils_source_objs, 'num_ils_parameters'))
            header_sect.Add_Child(Find_Keyword(ils_source_objs, 'Num_Ils_Wndepend'))
            try: nd_function_type
            except NameError:
                print "Couldn't find function_type.  Skipping..."
            else:
                header_sect.Add_Child(nd_function_type)
            header_sect.Add_Child(Find_Keyword(ils_source_objs, 'Ils_Cycle'))
            header_sect.Add_Child(Find_Keyword(ils_source_objs, 'Apo_M'))
            header_sect.Add_Child(Find_Keyword(ils_source_objs, 'Interpolation'))
            header_sect.Add_Child(L2_Input.Node('comment', '# resolution of fine grid/calculated grid'))
            Clean_Extra_Space(header_sect)
            Clean_Assignment_Spacing(header_sect)   

            fil_ils_out = "%s/%s" % (mod_dir, dst_ils_curr_name)
            Backup_Write(dst_ils_fileobj, fil_ils_out)

    # Move retrieval_mode from oco_l2.inp to aerosol files
    aerosol_ap_filenames=[]
    retrieval_mode = None
    for pd_sub_sect in param_def_sub_sects:
        if pd_sub_sect.leaf[0].lower() == 'aerosol':
            keyword_objs = pd_sub_sect.Get_All_Keyword_Objs()

            for key_obj in keyword_objs:
                if key_obj.leaf.lower() == 'a_priori':
                    aerosol_ap_filenames.append(key_obj.children[0].leaf)
                elif key_obj.leaf.lower() == 'first_guess':
                    aerosol_ap_filenames.append(key_obj.children[0].leaf)
                elif key_obj.leaf.lower() == 'retrieval_mode':
                    ret_mode = key_obj.children[0].leaf
                    if retrieval_mode != None:
                        if ret_mode != retrieval_mode:
                            raise ValueError('Inconsistent retrieval mode for AEROSOL')
                        else:
                            retrieval_mode = ret_mode

            pd_sub_sect.Delete_Keyword('retrieval_mode')
                
    for curr_aer_filename in aerosol_ap_filenames:
        #########################################################
        if os.path.exists(curr_aer_filename):
            print "Updating '" + curr_aer_filename + "'"
        else:
            print "Skipping missing: %s" % curr_aer_filename
            continue
        #########################################################

        dst_aer_fileobj = L2_Input.Input_File(curr_aer_filename)

        aer_source_objs = copy.copy(source_objs)
        aer_source_objs.append(dst_aer_fileobj)
            
        header_sect_srch = Find_Section([dst_ils_fileobj], 'HEADER')
        if header_sect_srch == None or len(header_sect_srch) == 0:
            raise IOError('Could not find section HEADER in file: %s' % curr_aer_filename)
        header_sect = header_sect_srch[0]
        header_sect.Set_Keyword('retrieval_mode', retrieval_mode)

        Clean_Extra_Space(header_sect)
        Clean_Assignment_Spacing(header_sect)   

        fil_aer_out = "%s/%s" % (mod_dir, curr_aer_filename)
        Backup_Write(dst_aer_fileobj, fil_aer_out)


    #########################################################
    print "Updating 'oco_l2.inp'"
    #########################################################
    
    # commit write of parameter section to file
    Clean_Assignment_Spacing(param_def_sect)   
    dst_inp_fileobj.rootNode.Add_Child(param_def_sect)
    Backup_Write(dst_inp_fileobj, "%s/%s" % (mod_dir, dst_inp_name))


    #########################################################
    print "Updating 'oco_l2.run'"
    #########################################################

    dst_run_fileobj = L2_Input.Input_File()

    # Create new control section
    control_sect = L2_Input.Section("CONTROL")
    control_sect.Add_Child(Find_Keyword(source_objs, 'input_file', 'oco_l2.run'))
    control_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    control_sect.Add_Child(Find_Keyword(source_objs, 'run_mode', 'RETRIEVAL'))
    control_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line
    
    control_sect.Add_Child(Find_Keyword(source_objs, 'jacobian_mode', 'ANALYTIC'))
    control_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line
    
    control_sect.Add_Child(Find_Keyword(source_objs, 'max_divergence'))
    control_sect.Add_Child(L2_Input.Node('comment', '# stop after this number of diverging steps'))

    control_sect.Add_Child(Find_Keyword(source_objs, 'max_iterations'))
    control_sect.Add_Child(L2_Input.Node('comment', '# stop after this many iterations'))

    control_sect.Add_Child(Find_Keyword(source_objs, 'max_chi2'))
    control_sect.Add_Child(L2_Input.Node('comment', '# fail convergences if chi2 > max_chi2'))

    control_sect.Add_Child(Find_Keyword(source_objs, 'lm_gamma'))
    control_sect.Add_Child(L2_Input.Node('comment', '# Levenberg-Marquardt gamma'))
    
    control_sect.Add_Child(Find_Keyword(source_objs, 'scale_convergence'))
    control_sect.Add_Child(L2_Input.Node('comment', '# converged if d_sigma_sq < scale_convergence * len_sv'))
    control_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line
    
    control_sect.Add_Child(Find_Keyword(source_objs, 'final_rad'))
    control_sect.Add_Child(Find_Keyword(source_objs, 'do_error'))
    # FAO -- fix following
    control_sect.Add_Child(Find_Keyword(source_objs, 'external_model_file', error_if_unfound=False))
    Clean_Extra_Space(control_sect)
    Clean_Assignment_Spacing(control_sect)
    dst_run_fileobj.rootNode.Add_Child(control_sect)

    # Create new algorithms section
    algs_sect = L2_Input.Section("ALGORITHMS")
    algs_sect.Add_Child(Find_Keyword(source_objs, 'polarization', 'true'))
    algs_sect.Add_Child(Find_Keyword(source_objs, 'num_diodes', '1016 1016 1016'))
    algs_sect.Add_Child(Find_Keyword(source_objs, 'streams', '8'))
    algs_sect.Add_Child(L2_Input.Node('comment', '# Number of streams that RADIANT will use'))
    #algs_sect.Add_Child(Find_Keyword(source_objs, 'interpolation'))
    #algs_sect.Add_Child(L2_Input.Node('comment', '# resolution of fine grid/calculated grid'))
    algs_sect.Add_Child(Find_Keyword(source_objs, 'points_sun'))
    algs_sect.Add_Child(Find_Keyword(source_objs, 'single_scatter_correction', 'true'))
    algs_sect.Add_Child(Find_Keyword(source_objs, 'delta_m_scaling', 'true'))
    algs_sect.Add_Child(Find_Keyword(source_objs, 'retrieval_algorithm', 'connor'))
    algs_sect.Add_Child(Find_Keyword(source_objs, 'h2o_correction', 'true'))
    algs_sect.Add_Child(Find_Keyword(source_objs, 'gas_interp_algorithm', 'sublayer'))
    Clean_Extra_Space(algs_sect)
    Clean_Assignment_Spacing(algs_sect)   
    dst_run_fileobj.rootNode.Add_Child(algs_sect)

    # Create new output section
    output_sect = L2_Input.Section("OUTPUT")
    output_sect.Add_Child(Find_Keyword(source_objs, 'log_file'))
    output_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    output_sect.Add_Child(L2_Input.Node('comment', '# Choices for verbosity are are DEBUG, INFO, WARNING, ERROR, FATAL.'))
    output_sect.Add_Child(Find_Keyword(source_objs, 'verbosity'))
    output_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    output_sect.Add_Child(Find_Keyword(source_objs, 'append'))
    output_sect.Add_Child(Find_Keyword(source_objs, 'control_flag'))
    output_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    output_sect.Add_Child(Find_Keyword(source_objs, 'output_path'))
    output_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    output_sect.Add_Child(L2_Input.Node('comment', '# Save jacobians and covariances in this dir for offline error'))
    output_sect.Add_Child(L2_Input.Node('comment', '# analysis.  If not present, no files will be saved.'))

    #$$$
    output_path = Find_Keyword(source_objs, 'output_path')[0].children[0].leaf
    nd_control_path = Find_Keyword(source_objs, 'control_path')
    nd_controlsub_path = Find_Keyword(source_objs, 'controlsub_path')
    nd_diagnostic_path = Find_Keyword(source_objs, 'diagnostic_path')

    control_path = output_path + nd_control_path[0].children[0].leaf
    controlsub_path = control_path + nd_controlsub_path[0].children[0].leaf
    diagnostic_path = output_path + nd_diagnostic_path[0].children[0].leaf

    nd_control_path[0].children[0].leaf = control_path
    nd_controlsub_path[0].children[0].leaf = controlsub_path
    nd_diagnostic_path[0].children[0].leaf = diagnostic_path

    output_sect.Add_Child(nd_control_path)
    output_sect.Add_Child(nd_controlsub_path)
    output_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    output_sect.Add_Child(L2_Input.Node('comment', '# All diagnostic files go here.  If not present, no files will be created.'))

    output_sect.Add_Child(nd_diagnostic_path)
    output_sect.Add_Child(L2_Input.Node('comment', ''))
    output_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    output_sect.Add_Child(Find_Keyword(source_objs, 'result_file'))
    output_sect.Add_Child(L2_Input.Node('comment', '# state vector output file name'))
    output_sect.Add_Child(Find_Keyword(source_objs, 'out_info_file'))
    output_sect.Add_Child(L2_Input.Node('comment', '# additional output file'))
    output_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    output_sect.Add_Child(Find_Keyword(source_objs, 'atmos_file'))
    output_sect.Add_Child(Find_Keyword(source_objs, 'summary_file', 'summary.dat'))
    output_sect.Add_Child(Find_Keyword(source_objs, 'control_file'))
    output_sect.Add_Child(L2_Input.Node('comment', '# in control_path'))
    output_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    output_sect.Add_Child(Find_Keyword(source_objs, 'output_each_iteration'))
    output_sect.Add_Child(L2_Input.Node('comment', '# save radiances & jacobians at each iteration'))
    output_sect.Add_Child(Find_Keyword(source_objs, 'high_res_spectra'))
    output_sect.Add_Child(L2_Input.Node('comment', '# create high-resolution spectra files'))
    output_sect.Add_Child(Find_Keyword(source_objs, 'save_jacobians'))
    output_sect.Add_Child(L2_Input.Node('comment', '# save jacobian files'))
    output_sect.Add_Child(Find_Keyword(source_objs, 'solar_transmittance', 'false'))
    output_sect.Add_Child(L2_Input.Node('comment', '\n')) # Empty line

    output_sect.Add_Child(Find_Keyword(source_objs, 'constraint_log_file', 'positive_constraint.log'))
    output_sect.Add_Child(Find_Keyword(source_objs, 'aggregator_dir', 'out/aggregator'))

    Clean_Extra_Space(output_sect)
    Clean_Assignment_Spacing(output_sect)   
    dst_run_fileobj.rootNode.Add_Child(output_sect)

    Backup_Write(dst_run_fileobj, "%s/%s" % (mod_dir, dst_run_name))


    #########################################################
    print "Updating 'oco_l2.win'"
    #########################################################
    
    dst_win_fileobj = L2_Input.Input_File()

    win_sects = Find_Section(source_objs, 'WINDOW')
    for curr_win_sect in win_sects:
        win_id = curr_win_sect.Get_Keyword_Value('id')
        if win_id.isdigit() and int(win_id) < len(new_window_id)+1:
            curr_win_sect.Set_Keyword_Value('id', new_window_id[int(win_id)-1])
        
        Clean_Extra_Space(curr_win_sect)
        Clean_Assignment_Spacing(curr_win_sect)
        dst_win_fileobj.rootNode.Add_Child(curr_win_sect)

    Backup_Write(dst_win_fileobj, "%s/%s" % (mod_dir, dst_win_name))

    ndlist = sound_info_sect.Get_All_Keyword_Objs()
    for nd in ndlist:
        if nd.leaf == 'soundinginfo_file':
            fil_soundinginfo = nd.children[0].leaf

    try: fil_soundinginfo
    except NameError:
        print "No sounding file specified in input; skip update."
    else:
        try:
            #########################################################
            print "Updating '" + fil_soundinginfo + "'"
            #########################################################

            # copy to data to buffer before we overwrite
            soundinfo_obj = L2_Input.Input_File(fil_soundinginfo)

            nd = Find_Keyword(source_objs, 'zero_azimuth', error_if_unfound=False)
            if len(nd) > 0:
                val_zero_azimuth = nd[0].children[0].leaf
            else:
                val_zero_azimuth = 'false'
                
            soundinfo_obj.Set_Keyword_Value('zero_azimuth', val_zero_azimuth)
            soundinfo_obj.Set_Keyword_Value('polarization_angle', '90.0D0')
            soundinfo_obj.Set_Keyword_Value('relative_velocity', '0.0D0')

            Clean_Extra_Space(soundinfo_obj.rootNode)
            Clean_Assignment_Spacing(soundinfo_obj.rootNode)

            Backup_Write(soundinfo_obj, fil_soundinginfo)
        except IOError:
            print "Couldn't find sounding file.  Skipping update."

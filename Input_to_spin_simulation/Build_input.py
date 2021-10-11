
import os
import shutil

def build_inputdat(source_path, new_path, params):
    old_file = open(source_path, "r")
    lines = old_file.read().split("\n")
    old_file.close()
    materials, widths, laser, closed = params
    
    if (len(materials) != len(widths)):
        raise ValueError("Materials and widths not matching")
        
    material_string = " ".join(materials)
    z_max = f"{int(sum(widths))}."
    N_regions = str(len(materials))
    # print(materials)
    # print(N_regions)
    laser  = f"{laser}."
    interfaces = []
    z = 0
    for w in widths[:-1]:
        z += w
        interfaces.append(str(z))
    interfaces = ". ".join(interfaces) + "."
    
    lines[12] = z_max
    lines[26] = N_regions
    lines[28] = interfaces
    lines[30] = material_string
    lines[38] = laser
    
#     for i, line in enumerate(lines):
#         print(f"{i}) " + line)
    
    
    new_file_string = "\n".join(lines)
    new_file = open(new_path,"w")
    new_file.write(new_file_string)
    new_file.close()

def build_submitsh(source_path, new_path, params):
    old_file = open(source_path, "r")
    lines = old_file.read().split("\n")
    old_file.close()
    materials, widths, laser, closed = params
    m = "".join(materials)
    w = "-".join([str(w) for w in widths])
    l = str(laser)
    if closed:
        c = "closed"
    else:
        c = "open"
    description = f"W{w}-L{l}"
    lines[1] += description
    
    new_file_string = "\n".join(lines)
    new_file = open(new_path,"w")
    new_file.write(new_file_string)
    new_file.close()

def generate_params(materials, widths, lasers, closed):
    """Generates a list of param tuples that define material, widths and laser"""
    param_list = []
    
    
    if len(materials) == 2:
        wlist1, wlist2 = widths
        for w1 in wlist1:
            for w2 in wlist2:
                for laser in lasers:
                    for c in closed:
                        param_list.append((materials, (w1,w2), laser, c))
    
    elif len(materials) == 3:
        wlist1, wlist2, wlist3 = widths
        for w1 in wlist1:
            for w2 in wlist2:
                for w3 in wlist3:
                    for laser in lasers:
                        for c in closed:
                            param_list.append((materials, (w1,w2,w3), laser, c))
    else:
        raise ValueError("2 or 3 material layers accepted")
    return param_list

def params_to_path(parent_path, params):
    materials, widths, laser, closed = params
    m = "".join(materials)
    w = "-".join([str(w) for w in widths])
    l = str(laser)
    if closed:
        c = "closed"
    else:
        c = "open"
    path = f"{parent_path}/{m}-{c}/W{w}-L{l}"
    return path

def build_one_folder(folder_path, source_path, params):
    
    try:
        shutil.rmtree(folder_path)
    except:
        pass
    
    os.makedirs(folder_path)
    materials_data_source_path = f"{source_path}/materialsdata_up_down.dat"
    materials_data_target_path = f"{folder_path}/materialsdata_up_down.dat"
    shutil.copyfile(materials_data_source_path, materials_data_target_path)
    
    closed = params[3]
    o, c = "open", "closed"
    reflect_source_path = f"{source_path}/reflect_ssv1_{c*int(closed)+o*int(not closed)}.dat"
    reflect_target_path = f"{folder_path}/reflect_ssv1.dat"
    shutil.copyfile(reflect_source_path, reflect_target_path)
    
    input_source_path = f"{source_path}/input.dat"
    input_target_path = f"{folder_path}/input.dat"
    build_inputdat(input_source_path, input_target_path, params)
    
    submit_source_path = f"{source_path}/submit.sh"
    submit_target_path = f"{folder_path}/submit.sh"
    build_submitsh(submit_source_path, submit_target_path, params)

def build_all_folders(param_combos, parent_path, source_path):
    materials, widths, laser, closed = param_combos
    params_list = generate_params(materials, widths, laser, closed)
    for params in params_list:
        input_path = params_to_path(parent_path, params)
        build_one_folder(input_path, source_path, params)



if __name__ == "__main__":

    materials= ["Fu", "Pt"]
    w1 = list(range(2,7))
    w2 = list(range(2,7))
    widths = [w1, w2]
    laser = [30,60,100]
    closed = [True]

    param_combos = (materials, widths, laser, closed)

    parent_path = "/home/hjaleta/spin_project/define_input/generated_input"
    source_path = "/home/hjaleta/spin_project/define_input/input_source_files"

    build_all_folders(param_combos, parent_path, source_path)
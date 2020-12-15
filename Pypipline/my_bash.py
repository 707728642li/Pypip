import os
import sys
import re
import argparse
import platform
from collections import OrderedDict,defaultdict

# remove log file
rm_test_list = [i for i in os.listdir() if "log" in i]
for each in rm_test_list:
    os.remove(each)
    
# get command
on_server = True if platform.platform().lower().startswith("linux") else False
if on_server:
    parser = argparse.ArgumentParser(description="Pipysh = Python make shell pipline!")
    parser.add_argument("file",help="config file")
    parser.add_argument("-a","--action",
                        choices=["run","detail","basic"],default="basic"
                        ,help=(f"basic - show basic running infomation (default);"
                               f"detail - show detail running infomation;"
                               f"run - run bash file."
                                ))
    parser.add_argument("-t","--test"
                        ,default=":"
                        ,help="Program name:sample name")


    args = parser.parse_args()

    configure_file = args.file
    action = args.action
    t_items = args.test.split(":")
else:
    class args:
        test=":"
    configure_file = "bash_cmd2.txt"
    action = "basic"
    
class bio_pipline:
    total_log = "total_info.log.txt"
    
    @staticmethod
    def check_config(content):
        mistake_list = []
        for line in content:
            comma_n = line.split("->",1)[0].count(",")
            arrow_n = line.split(",",1)[1].count("->")
            for n,marker in zip((comma_n,arrow_n)
                                ,((",","before","->"),("->","after",","))):
                if n!=1:
                    mistake_list.append((n,marker,line))
#         print(mistake_list)
        if mistake_list:
            print("Input errors was found in these lines:")
            for each in mistake_list:
                print(f"{each[0]} '{each[1][0]}' was found {each[1][1]} '{each[1][2]}' in line: {each[2]}")
            raise ValueError("Input error!/(ㄒoㄒ)/~~")
            
    @staticmethod
    def get_dict(configure_file):
        """
        get infomation from configure file
        1.2,p -> bwa mem genome.fa -a xxx_1.fq xxx_2.fq > xxx.sam
        1.2: program number . subgroup number : p_number
        p: commad type : c_type, including:
            i(xxx):input samples, separated by ";" and "\n", not required
              xxx: default sample label, xxx will be replaced by sample names
            t: subgroup number, 10 samples / 3 groups means: (4, 3, 3)
            p: process shell command
            
        """
        import os
        import re
        from collections import OrderedDict#,defaultdict
        config_dict = OrderedDict()
        #read configure file
        all_content = open(configure_file).readlines()
        all_content = [i for i in all_content if i.strip()]
        
        bio_pipline.Stop_if_error = " ; "
        error_treat = [i for i in all_content 
                       if i.startswith("#") and "stop_if_error" in i.lower()]
        
        if error_treat and error_treat[0].strip().split()[-1].lower() in ["true","t"]:
                bio_pipline.Stop_if_error = " && "
             
        all_content = [i for i in all_content if not i.startswith("#")]
        sub_p_dict = defaultdict(set)
        program_seq = []
        
        bio_pipline.check_config(all_content)

        for line in all_content:
            # program number
            p_number = line.split(",")[0].strip(".")
            
            if "." in p_number and p_number.split(".")[1].strip():
                p_mother_name = p_number.split(".")[0]
                sub_p_dict[p_mother_name].add(p_number)
            else:
                p_mother_name = p_number
            if not (p_mother_name in program_seq):
                program_seq.append(p_mother_name)
                
            config_dict.setdefault(p_number,
                                   {"sample_label":"xxx",
                                    "t":"1",
                                    "i":set()})

            # detect i p t
            c_type = line.split(",")[1].strip()[0]
            if c_type=="i":
                # detect sample label, default:xxx
                try:
                    sample_label = re.findall(r"i\(\s*(.*?)\s*\)",line.split("->")[0])[0]
                    if sample_label:
                        config_dict[p_number]["sample_label"] = sample_label
                except:
                    pass
                #if sample_label:
                #    config_dict[p_number]["sample_label"] = sample_label
            
                # add input samples
                content = line.split("->",1)[1]
                if "from-file=" in content.lower():
                    in_file_name = content.split("=")[-1].strip()
                    content = open(in_file_name).readlines()
                    content = [i.split(";") for i in content]
                    content = [j.strip() for i in content for j in i if j.strip()]
#                     print(content)
                else:
                    content = content.strip().strip(";").split(";")
                    # remove all space in two sides for sample_name
                    content = [i.strip() for i in content if i]
                config_dict[p_number][c_type].update(content)

            if c_type=="t":
                threads = line.split("->")[1].strip()
                config_dict[p_number]["t"] = threads

            if c_type == "p":
                # split parameters
                content = line.split("->",1)[1].strip(";").split(";")
                # remove all space in two sides
                content = [i.strip() for i in content]
                config_dict[p_number].setdefault(c_type,[]).extend(content)
        no_label = []
        for p in config_dict.keys():
            config_dict[p]["i"] = sorted(list(config_dict[p]["i"]))
            config_dict[p]["i"] = [i for i in config_dict[p]["i"] if i]
            config_dict[p]["i"] = ["just.shell.command.without.sample.name"] if not config_dict[p]["i"] else config_dict[p]["i"]
            config_dict[p]["t"] = config_dict[p].get("t","1") 
            
            if config_dict[p]["i"] != ["just.shell.command.without.sample.name"]:
                sample_label = config_dict[p]["sample_label"]
                no_label.extend([(p,sample_label,i) for i in config_dict[p]["p"] if sample_label not in i])
        f_warning = lambda x: f"Program[{x[0]}](label-{x[1]}): {x[2]}"
        warning_info = ["No labels were found in these commands:"]+[f_warning(i) for i in no_label]
#         print(warning_info)
        if no_label:
            bio_pipline.print_head(warning_info,l=75) 

        bio_pipline.config_dict = config_dict # 包含 所有信息
        bio_pipline.program_subdict = sub_p_dict # 包含 项目.的子项目信息
        bio_pipline.program_seq = program_seq # 包含 运行顺序信息
        
    #============================= Prepare samples =============================#
    
    # split samples into different groups
    @staticmethod
    def split_sample(sample_list,n=1):
        #print(sample_list)
        sample_list = list(sample_list)[::-1]
        out_dict = {}
        while sample_list:
            for i in range(1,n+1):
                #print(i)
                if sample_list:
                    out_dict.setdefault(i,[]).append(sample_list.pop())
                if not sample_list:
                    break
        return out_dict

    # get commands for group
    @staticmethod
    def get_cmd_per_group(label,sample_list,cmd,):
        out_cmd_list = []
        sample_list = list(sample_list)
        out_dict = {}
        #error mark: " ; " ignore=True else " && "
        e_mark = bio_pipline.Stop_if_error
        for each in sample_list:
            if each == "just.shell.command.without.sample.name":
                out_dict[each] = e_mark.join(cmd)
            else:
                c_each = f"({each})" if show_only else each
                # if run or test, no '()'
                out_dict[each] =  e_mark.join([i.replace(label,c_each) for i in cmd])

        return out_dict #{sample_1:cmd_1,sample_2:cmd_2}
    
    # use above method to get commad_dict for each p 
    @staticmethod
    def get_cmd_per_program(p_dict):
        sample = tuple(p_dict["i"])
        try:
            threads = int(p_dict["t"]) 
        except:
            threads = len(sample) if p_dict["t"]=="!" else 1
            
        sample_dict = bio_pipline.split_sample(sample,threads)
        label = p_dict["sample_label"]
        cmd = p_dict["p"]
        run_dict = {}
        for sub_num,sub_list in sample_dict.items():
            run_dict[sub_num] = bio_pipline.get_cmd_per_group(label,sub_list,cmd)
#         print(run_dict)

        return run_dict #{tread_1:{sample_1:cmd_1,sample_2:cmd_2}}
        
    #============================= Run methods =============================#

    @staticmethod
    def run_popen(p_name,subgroup_n,cmd_dict):
        """
        run popen for each subgoup under each program
        call by Pool() in run_cmd_per_program()
        p_name: program name
        subgroup_n: subgroup_n number, or subgroup number, from 1
        cmd_dict: cmd dict for each sample in this group   
        """
        import subprocess
        import os
        import time
        from datetime import datetime
        if show_only:
            for sample,cmd in cmd_dict.items():
                sample = "Shell_command" if sample =="just.shell.command.without.sample.name" else sample
                print(f"Sample: {sample}:  {cmd}")            
        else:
            for sample,cmd in cmd_dict.items():
                if sample!="just.shell.command.without.sample.name":
                    p_log = (f"Program-{p_name}_group-{subgroup_n}_sample-{sample.replace(' ','_')}.log.txt")
                else:
                    sample = "Shell_command"
                    p_log = (f"Program-{p_name}.log.txt")
                print(sample)
                if os.path.exists(p_log):
                    os.remove(p_log)
                with open(p_log,"w") as log_file:
                    start_time = datetime.now()
                    start_info = (f"Sample: {sample}"
                        f"\nStart time: {str(start_time)}\n{cmd}\n")
                    log_file.write(start_info)
                    print(f"Start to run: {p_log.rstrip('.log.txt')}")
                    
                with open(p_log,"a") as log_file:
                    # running...
                    subjob = subprocess.Popen(cmd,shell=True,
                                              bufsize=1,universal_newlines=True,
                                             stdout=log_file,stderr=log_file)
                    subjob.wait()
                    
                    end_time = datetime.now()
                    cost_time = round((end_time - start_time).total_seconds()/3600,2)
                    log_file.write(f"\nFinished!\n{str(end_time)}\nTime costs: {cost_time} hours\n")
                if subjob.returncode ==0:
                    return_info = f"Successfully!(*^_^*) Time costs: {cost_time} hours" 
                else:
                    return_info = f"Failed!/(ㄒoㄒ)/~~ Time costs: {cost_time} hours. Please check log file: {p_log}" 
                print((f"Program {p_name} - {sample}: {return_info}"))
                with open(bio_pipline.total_log,"a") as total_log:
                    total_log.write(f"Program {p_name} - {sample}: {return_info}\n\n")

    #============================= Run samples =============================#

    @staticmethod    
    def run_cmd_per_program(p_name,run_dict,):
        """
        use multiprocessing.Pool(run_popen()) to run each programs
        run_dict = get_cmd_per_program()
        
        """
        if show_only:
            for subgroup_n,cmd in run_dict.items():
                n_sample = len(cmd.keys()) if list(cmd.keys())!=["just.shell.command.without.sample.name"] else 0
                print(f"Sub group {subgroup_n}: with {n_sample} samples")
                bio_pipline.run_popen(p_name,subgroup_n,cmd
                                      ,)
        else:
            from multiprocessing import Pool
            n = len(run_dict.keys())
            pool = Pool(n)
            for subgroup_n,cmd in run_dict.items():
                pool.apply_async(
                    bio_pipline.run_popen,
                    (p_name,subgroup_n,cmd,))
            pool.close()
            pool.join()

    @staticmethod
    def run_single_work(p_name,):
        """
        run for single programs, programs without "."
        """
        job_dict = bio_pipline.config_dict[p_name]
        job_cmd_dict = bio_pipline.get_cmd_per_program(job_dict)
        if show_only:
            print(f"\nRun Program [{p_name}]")
        bio_pipline.run_cmd_per_program(p_name,job_cmd_dict
                                        ,)
    @staticmethod
    def run_sub_works(sub_work_set,):
            
        """
        run for multiple programs, programs with "."
        call run_single_work()
        
        """
        from multiprocessing import Pool
        sub_work_list = sorted(list(sub_work_set))
        sub_n = len(sub_work_list)
        pool = Pool(sub_n)
        print(f"\nRun Subprogram {' '.join([f'[{i}]' for i in sub_work_list])} in same time\n")
        for sub in sub_work_list:
            if show_only:
                bio_pipline.run_single_work(sub,)
            else:
                pool.apply_async(bio_pipline.run_single_work,(sub,))
                pool.close()
                pool.join()
        
    @staticmethod
    def run_main(configure_file,):
        bio_pipline.get_dict(configure_file)
        for each in bio_pipline.program_seq:
            if each not in bio_pipline.program_subdict.keys():
                bio_pipline.run_single_work(each,)
            else:
                bio_pipline.run_sub_works(bio_pipline.program_subdict[each],)
        
    @staticmethod
    def run_test(configure_file,):
        bio_pipline.get_dict(configure_file)
        t_program = t_items[0]
        if len(t_items) == 2:
            t_sample = t_items[1]
        else:
            # Just use the first sample
            t_sample = bio_pipline.config_dict[t_program]['i'][0]
            
        if t_program not in bio_pipline.config_dict.keys():
            raise ValueError(f"\nProgram name {t_program} was not found in Program names:"
                            f"\n{' '.join(bio_pipline.config_dict.keys())}")
        if t_sample not in bio_pipline.config_dict[t_program]["i"]:
            raise ValueError(f"\nSample name {t_sample} was not found in Program names:"
                            f"\n{' '.join(bio_pipline.config_dict[t_program]['i'])}")
        bio_pipline.config_dict[t_program]["i"] = [t_sample]
        bio_pipline.print_head(f"Command test for Program: [{t_program}] and sample: {t_sample}")
        bio_pipline.run_single_work(t_program,)
        
    #============================= Show infomations =============================#
    @staticmethod
    def print_head(s,l=60):
        u_d = "‖=" + f"".center(l,"=") + "=‖"
#         u_d = " =" + f"".center(l,"=") + " ‖"
        c_str = "‖ " + f"{s}".center(l," ") + " ‖"
        if type(s)==str:
            out_str = u_d+"\n"+ c_str + "\n" + u_d + "\n"
            print(out_str)
        else:
            out_content = "\n".join(["‖ " + f"{i}".ljust(l," ") + " ‖" for i in s])
            out_str = u_d+"\n"+ out_content + "\n" + u_d
            print(out_str)
    @staticmethod
    def show_single_program(p_name):
        s_dict = bio_pipline.config_dict
        print("\n"+f" Program [{p_name}] ".center(80,"="))
#         bio_pipline.print_head(f"Program [{p_name}]",l=80)
        empty_sample = True if s_dict[p_name]["i"] == ['just.shell.command.without.sample.name'] else False
        n_sample = len(s_dict[p_name]['i']) if not empty_sample else 0
        print(f"{n_sample} samples in {s_dict[p_name]['t']} groups")
        for i,j in bio_pipline.split_sample(s_dict[p_name]['i'],int(s_dict[p_name]['t'])).items():
            print(f"Subgroup {i} with {len(j) if not empty_sample else 0} samples: {' ; '.join(j)}")
        temp_label = s_dict[p_name]['sample_label'] if not empty_sample else ""
        print(f"\nCommand↓ Sample Label: {temp_label}")
    #     print(s_dict[p_name]["i"])
        for each in s_dict[p_name]["p"]:
            if not empty_sample:
                each = each.replace(temp_label,f"({temp_label})")
            if each.count(temp_label)==0 and s_dict[p_name]["i"] != ['just.shell.command.without.sample.name']:
                each = each + f"    ★ no '{temp_label}' find in this command ★"
            print(each)
        print("="*80)
        
    @staticmethod
    def show_basic():
        bio_pipline.get_dict(configure_file)
        stop = "True" if bio_pipline.Stop_if_error == " && " else "False"
        print(f"Stop program if error: {stop}")
        for each in bio_pipline.program_seq:
            if each not in bio_pipline.program_subdict.keys():
                bio_pipline.show_single_program(each)
    #             bio_pipline.run_single_work(each,)
            else:
                print("\n"+f" Run Subprogram {' & '.join(bio_pipline.program_subdict[each])} in same time ".center(80,"*"))
                sub_list = sorted(list(bio_pipline.program_subdict[each]))
                for sub in sub_list:
                    bio_pipline.show_single_program(sub)
                print("\n"+"".center(80,"*"))


    #============================= Run! =============================#
if __name__ == "__main__":
    if args.test != ":":
        show_only = False
        bio_pipline.run_test(configure_file=configure_file)
    else:
        if action == "basic":
            show_only = True
            bio_pipline.show_basic()
        elif action in {"run","detail"}:
            show_only = True if action=="detail" else False          
            bio_pipline.run_main(configure_file=configure_file,)
        else:
            print("Input Error")


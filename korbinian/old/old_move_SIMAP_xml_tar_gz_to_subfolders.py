import os
import korbinian
from time import strftime, gmtime

#get the time, including seconds
the_time = strftime("%Y-%m-%d_%H_%M_%S")
#logfile = r"D:\Databases\simap\test_distribute\move_simap_logfile.txt"
logfile = "/nas/teeselab/mark/files/main/logfiles/{}_move_simap_rimma_to_db_logfile.txt".format(the_time)
logging = korbinian.mtutils.setup_error_logging(logfile)

# source_dir = r"D:\Databases\simap\test_distribute"
# source_dir = r"D:\Databases\simap\Archaea"
# source_dir = r"D:\Databases\simap\Bacteria"
source_dir = "/nas/teeselab/students/rimma/databases/simap/00_test"
target_dir = "/nas/teeselab/db/simap"

def absoluteFilePaths(directory):
    # http://stackoverflow.com/questions/9816816/relative-and-absolute-paths-of-all-files
    for dirpath, _, filenames in os.walk(directory):
        for f in filenames:
            yield os.path.abspath(os.path.join(dirpath, f))

list_all_filepaths = list(absoluteFilePaths(source_dir))

logging.info("\nsource_dir:\n{s}\n{t}\n".format(s=source_dir, t=target_dir))

for filepath in list_all_filepaths:
    filedir = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    last_13 = filename[-13:]
    if last_13 == "_SIMAP.tar.gz":
        first_2 = filename[:2]
        target_subdir = os.path.join(target_dir, first_2)
        if not os.path.exists(target_subdir):
            #print("folder does not exist:\n{}".format(target_subdir))
            os.makedirs(target_subdir)
        else:
            pass
            #print("folder exists:\n{}".format(target_subdir))
        path_orig = os.path.join(filedir, filename)
        logging.info("path_orig  : {}".format(str(path_orig)))
        path_target = os.path.join(target_dir, first_2, filename)
        logging.info("path_target: {}".format(str(path_target)))
        if not os.path.isfile(path_target):
            # move file
            os.rename(path_orig, path_target)
            logging.info("{} moved".format(filename[:6]))
        else:
            logging.info("{} already in target folder".format(filename[:6]))
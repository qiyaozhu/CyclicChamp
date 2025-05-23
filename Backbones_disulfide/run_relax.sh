source /mnt/home/vmulligan/load_my_modules.sh
source /mnt/home/vmulligan/masala_workingcopy/set_up_masala.sh

/mnt/home/vmulligan/rosetta_git_workingcopy/Rosetta/main/source/bin/rosetta_scripts.masalampiserialization.linuxgccrelease -masala_plugins ~vmulligan/masala_workingcopy/standard_masala_plugins/ -total_threads 4 -masala_total_threads 4 @relax.flags

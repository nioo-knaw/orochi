#!/usr/bin/python
assembly_mode = str(input("What type of assembly would you like to use? (1 - all; 2 - per treatment; 3 - per sample) "))
if assembly_mode == "1":
    import Configmaker_all_1
if assembly_mode == "2":
    import Configmaker_treatments_2
if assembly_mode == "3":
    import Configmaker_samples_3

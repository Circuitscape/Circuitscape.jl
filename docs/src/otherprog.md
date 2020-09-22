<<<<<<< HEAD
# Calling Circuitscape from other programs

Circuitscape can be invoked from external programs and scripts (e.g., Python and R) to do computations on rasters and networks and return results.

More soon.
=======
# 7\. Calling Circuitscape from other programs

Circuitscape can be invoked from external programs and scripts (e.g., Python and R) to do computations on rasters and networks and return results. It reads user settings from a configuration (.ini) file that can be either created manually or saved from the user interface (under File >> Save settings).

If the external program can interface with Python packages, then the Circuitscape methods can be called directly when Circuitscape is installed as a Python package. For example:

    from circuitscape import Compute

    cs = Compute('configuration.ini', 'Screen')
    result = cs.compute()

If such direct interfacing with Python packages is not possible, Circuitscape can be invoked as an application with the following command:

    python csrun.py [configuration.ini]

Alternatively, the cs_run.exe executable can be invoked from external programs and scripts on Windows. Linkage Mapper calls Circuitscape this way from the Pinchpoint and Centrality modules, as does the Circuitscape for ArcGIS toolbox.

In the cases above, results are written out to files which are then read back by the external program.
>>>>>>> b4c223922c4b8ceaa698e62737a96d1eba542d29

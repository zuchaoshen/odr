## Resubmission

This is for updating the package. In this version, we added the ant colony 
optimization algorithm to identify optimal sample allocations for od.2m and 
od.3m functions. We also added functions to identify sample allocations for 
expeirments investigating mediation and moderation effects.

Across operational systems we tested, there are two notes that 
can be ignored (detailed below).

Please help process this package at your convenience! Thanks a lot! 

## Test environments
* local Windows 11, R 4.3.1
* r-hub
* win-builder (devel and release)

## R CMD check results
0 errors √ | 0 warnings √ | 0 notes √ (local Windows 11)
0 errors √ | 0 warnings √ | 2 notes ✖  (R-hub)
0 errors √ | 0 warnings √ | 0 note √  (win-builder )

## Reverse dependencies

There are no issues on reverse dependencies.

---
checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ''NULL''
  Found the following files/directories:
    'lastMiKTeXException'
    
 Response: This may be associated with the R-hub testing environment (Windows Server 2022, R-devel, 64 bit).

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found

Response: This may be associated with the R-hub testing environments (Fedora Linux, R-devel, clang, gfortran; Ubuntu Linux 20.04.1 LTS, R-release, GCC).
 

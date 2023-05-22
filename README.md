# monoRTM

---
**Contents**

1. [Cloning](#cloning)
2. [Directory Structure](#structure)
3. [Building monoRTM](#build)
4. [Running monoRTM](#run)

If any build or run issues occur, please [create an issue](https://github.com/AER-RC/monoRTM/issues) or contact the [AER-RC Group](https://github.com/AER-RC).

# Cloning the Latest Release <a name="cloning"></a>

Assuming the output directory should be `monoRTM`:

```
git clone --recursive git@github.com:AER-RC/monoRTM.git
```

`--recursive` is important, because this repository is linked with our [common FORTRAN modules repository](https://github.com/AER-RC/aer_rt_utils) that are required in the model builds. If this keyword is forgotten, one can do:

```
git submodule init
git submodule update
```

in the `monoRTM` directory.

Currently, the latest release is monoRTM v5.6, and it is recommended that this be the version that users clone and checkout (rather than the `master` branch). To do this, one needs to simply checkout the `v5.6` tag:

```
git checkout tags/v5.6
```

No releases before v5.6 are available via GitHub, but they can be requested by emailing <aer_monortm@aer.com>. For information on previous releases, please visit the [What's New Wiki page](https://github.com/AER-RC/monoRTM/wiki/What's-New).

Instead of cloning, users can also download a monoRTM [tarball](https://github.com/AER-RC/monoRTM/releases/tag/v5.6) and unpack it:

```
tar xvf monortm_v5.6.tar.gz
```

# Directory Structure <a name="structure"></a>

The `MONORTM` directory contains several sub-directories described briefly below:

| Name | Description |
| :--- | :--- |
| `README` | This file |
| `build` | Contains makefiles for MonoRTM for different platforms |
| `run` | Contains files required to run MonoRTM.<ul><li>`run_monortm_example` is a script to run some examples</li><li>`TAPE3_spectral_lines.dat.0_55.v5.0_fast` contains the spectral line information<ul><li>Note that if you wish to run a faster calculation with fewer microwave spectral lines or for a spectral range outside the wavenumber range of the spectral line file supplied in the examples, you can create a line file using MonoLNFL, which iscontained in this release package.</li></ul></li><li>`run/in` contains sample input files:<ul><li>`_dn`: denotes sample input files for downwelling radiance calculations</li><li>`_up`: denotes sample input files for upwelling radiance calculations</li><li>`MONORTM.IN`: input compatible with lblrtm (ex-`TAPE5`), `IATM=1`. IDL code to generate profiles in `TAPE5` format from ARM netCDF radiosonde files is available in the `idl` sub-directory</li><li>`MONORTM_PROF.IN` : contains layer data, `IATM=0` (this is a copy of `TAPE7` which is generated by a run where `IATM=1`</li></ul></li><li>`run/out` will hold the output from MonoRTM:<ul><li>`MONORTM.OUT`</li><li>`MONORTM.LOG` (if LBLATM is ON) - `TAPE6` file</li><li>`TAPE7` (if LBLATM is ON)</li></ul>
| `src` | contains all source files needed by MonoRTM. |
| `doc` | detailed instructions manual about MonoRTM, in ASCII format. |
| `idl` | tool for creating monortm input from ARM sonde files. |

# Building monoRTM <a name="build"></a>

To start, descend into the `build` directory:

```
cd build
make -f make_monortm $TARGET
```

The `TARGET` environment variable depends on the user's operating system, compiler, and desired precision. Available targets are:

| Target | Description | Compiler |
| :--- | :--- | :--- |
| `aixIBMsgl` | IBM/AIX OS using IBM fortran,single precision| `xlf90` |
| `linuxPGIsgl` | Linux OS using PGI fortran,single precision |  `pgf90` |
| `linuxGNUsgl` | Linux OS using GNU fortran,single precision | `gfortran` |
| `linuxG95sgl` | Linux OS using G95 fortran,single precision | `g95` |
| `inuxINTELsgl` | Linux OS using Intel fortran,single precision | `ifort` |
| `mingwGNUsgl` | Windows unix shell environment using gfortran,single precision | `gfortran` |
| `osxABSOFTsgl` | Mac OS_X using Absoft Pro fortran,singleprecision | `f90` |
| `osxGNUsgl` | Mac OS_X using GNU fortran,singleprecision | `gfortran` |
| `osxIBMsgl` | Mac OS_X using IBM XL fortran,singleprecision | `xlf90` |
| `osxINTELsgl` | Mac OS_X using Intel fortran,single precision | `ifort` |
| `sunSUNsgl` | Sun/Solaris OS using Sun fortran,single precision | `sunf90` |
| `sgiMIPSsgl` | SGI/IRIX64 OS using MIPS fortran,single precision | `f90` |
| `aixIBMdbl` | IBM/AIX OS using IBM fortran,double precision| `xlf90` |
| `linuxPGIdbl` | Linux OS using PGI fortran,double precision |  `pgf90` |
| `linuxGNUdbl` | Linux OS using GNU fortran,double precision | `gfortran` |
| `linuxG95dbl` | Linux OS using G95 fortran,double precision | `g95` |
| `inuxINTELdbl` | Linux OS using Intel fortran,double precision | `ifort` |
| `mingwGNUdbl` | Windows unix shell environment using gfortran,double precision | `gfortran` |
| `osxABSOFTdbl` | Mac OS_X using Absoft Pro fortran,doubleprecision | `f90` |
| `osxGNUdbl` | Mac OS_X using GNU fortran,doubleprecision | `gfortran` |
| `osxIBMdbl` | Mac OS_X using IBM XL fortran,doubleprecision | `xlf90` |
| `osxINTELdbl` | Mac OS_X using Intel fortran,double precision | `ifort` |
| `sunSUNdbl` | Sun/Solaris OS using Sun fortran,double precision | `sunf90` |
| `sgiMIPSdbl` | SGI/IRIX64 OS using MIPS fortran,double precision | `f90` |

# Running monoRTM <a name="run"></a>

To generate example output:

```
cd run
run_monortm_examples
```
Note that you may have to modify the executable name in `run_monortm_examples` to match the name of the executable that you have created.

Example inputs and outputs can be found in the [example tarball](https://github.com/AER-RC/monoRTM/files/6804136/monortm_examples.tar.gz).
```

```

# Tailoring monoRTM to Specific Needs

MONORTM is a driver program that calls the core module called the Monochromatic Optical Depth Model (MODM).  The inputs to MODM could be modified inside `monortm.f90` directly or through the input file `MONORTM.IN` (same type of format as LBLRTM's `TAPE5`), see instructions for more details.  MonoRTM is a forward model. MonoRTM is designed to be very flexible. We can either use it as a black box and control everything from the `MONORTM.IN` input file, or one can modify the code itself and recompile it. In the latter, it is structured in such a way that the changes should always be done in `monortm.f90` (the driver program). The other auxillary files (`monortm_sub.f90`, `modm.f90`, `lblatm.f90`, `lblrtm_sub.f`, `tips_2003.f` and `isotope.dat`) should normally not be touched, except in rare situations. In case there is an update in the continuum calculations a new `contnm.f90` module will be generated and sent to the users (or made available
on the WEB/ftp site).

The core of MonoRTM is the computation of the optical depths. It is designed as a subroutine for flexibility and could be easily plugged in a different radiative transfer model if needed.

# Outputs of MonoRTM:

The records stored in `MONORTM.OUT` are the following:

| Record | Description |
| :---: | :--- |
| `NPR`| Profile index used (not necessarily in order) |
| `FREQ` | Frequency in GHz (or wavenumbers, for wavenumbers greater than 100 cm<sup>-1</sup>) |
| `BT(I)`	| Brightness temperature in Kelvin |
| `TMR(I)` | Mean radiating temperature (K) |
| `RAD(I)` | Radiance (W/(cm^2 ster cm^-1)) |
| `TRTOT(I)` | Total transmittance (no unit: between 0 and 1) |
| `WVCOLMN` | Integrated water vapor amount along the optical path in cm |
| `CLWCOLMN` | Integrated cloud liquid water along the optical path in mm |
| `TMPSFC` | Surface/target temperature in K |
| `EMISS(I)` | Surface/target emissivity (no unit, between 0 and 1) |
| `REFLC(I)` | Surface/target Reflectivity (no unit , between 0 and 1) |
| `ANGLE`	| Angle in degrees |
| `OTOT` | Total column-integrated optical depth due to all species |
| `OTOT_*` | Total column-integrated optical depth by molecules with line data |
| `ODXSEC` | Total column-integrated optical depth due to all cross-section molecules |

# Contact

Any comments or questions should be forwarded to Karen Cady-Pereira (aer_monortm@aer.com)

AER , Radiation and Climate Group

131 Hartwell Avenue

Lexington, MA 02421 USA

Tel: 1 781 761 2216

Fax: 1 781 761 2299

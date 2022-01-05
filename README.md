# AVO Attributes
Program still under development.

Calculates and display Intercept, Gradient, Pseudo Poisson and Fluid Factor attributes given CRP gathers and a Velocity field file (ASCII). More attributes to be added in the future.

### Table of Contents
* [Installation](https://github.com/leorodve/AVO#installation)
  * [Linux](https://github.com/leorodve/AVO#linux)
    * [Install](https://github.com/leorodve/AVO#install)
* [Examples](https://github.com/leorodve/AVO#examples)
* [License & Contributing](https://github.com/leorodve/AVO#license-and-contributing)
* [Authors](https://github.com/leorodve/AVO#Authors)
* [References](https://github.com/leorodve/AVO#references)

## Installation
### Linux
### Install
You can download the latest version of the source code from [GitHub](https://github.com/leorodve/AVO) or if you have ```git``` installed, you can use the following command:
```bash
git clone https://github.com/leorodve/AVO.git
```
Change your directory to the newly created containing the source code and run:
```bash
python functions_AVO.py
```
## Examples
This program was tested using a seismic section containing multiple reflectors, one of which experience a decrease of its amplitude with the offset.

<img src="https://i.imgur.com/L0prznN.png">

The maximum angle of incidence is close to 20°, making it acceptable to use Shuey's two term AVO equation.

<img src="https://i.imgur.com/LqQ9Ige.png">

After uploading both files (CRP gathers and velocity field) the program calculates and displays the Intercept (Rp), Gradient (G), Pseudo Poisson (Δσ) and Fluid Factor(ΔF) using linear regression and Shuey's two term AVO equation.

<img src="https://i.imgur.com/JH7PMcv.png">

If we zoom in on the reflectors around 1300 ms we get:

<img src="https://i.imgur.com/ufqes8U.png">

While the intercept section shows the sandstone (red) with shale interbeddings (pale green) within the reservoir unit, the gradient section indicates variations in fluid saturation (indicated by the different tones of orange) within the reservoir rocks. The orange color in the fluid-factor section represents the gas-saturated reservoir sandstone.

## License and contributing
Licensed under the [GPL-3.0](http://www.gnu.org/licenses/gpl-3.0.html) License.

All kinds of contributions, including code, bug reports and issues are welcomed.

## Authors
Made by Leonardo Rodriguez.

## References
Yilmaz, O (2001). Seismic Data Analysis: Processing, Inversion and Interpretation of Seismic Data: 2nd (second) Edition. SEG

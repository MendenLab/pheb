# **P**redictive **H**TS **E**pigenetic **B**iomarkers (PHEB)

Environment
--- 
The analysis can be run in a local virtual [docker](https://docs.docker.com/) environment. For that, first pull the image by typing:
```
docker pull aljoshoh/pheb_r3.5.3
```
Finally, run the image by:
```
docker run -it -p 3838:3838 aljoshoh/pheb_r3.5.3 /bin/bash
```
You should now be in the bash of the pheb docker environment (if you are not sure, check for the existence of `/.dockerenv`). Now type:
```
start_rstudio 3838
```
Now you can navigate to `127.0.0.1:3838` in a browser for running a virtual rstudio session. If you want to build your own container, enter the env directory (`cd environment`) and run `docker build -f Dockerfile -t "aljoshoh/pheb_r3.5.3" .` Then, the above command can be used to run the generated docker image.

>If you are working on **Windows 10** you will have to convert the line endings of the two bash scripts. Herefore natigate to 
>```
>cd ..
>/usr/bin/
>```  
>and convert the line endings of your 2 bash scripts. 
>```
>sed -i -e 's/\r$//' start_rstudio.sh
>sed -i -e 's/\r$//' rstudio_auth.sh
>```
 
Now you can connect to a rstudio rsever session for costumization and running scripts manually. This will show you a link and password you can use to login on a webbrowser. You can nagivate to the link on a browser and login with username ***rstudio*** and the shown password.


Alternatively, you can use [conda](https://docs.conda.io/en/latest/) for installing its requirements. In the desired location, type:
```
git clone https://github.com/aljoshoh/pheb.git
cd pheb/
condo activate
./environment/make_env.sh
``` 
This will install a virtual conda environment called `pheb_r3.5.3`. It also contains an installation of `Rstudio`.

Analysis
---
The chunks provided in `workflow.Rmd` can be used to run the different substeps of the analysis.


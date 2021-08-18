# GroundContamination


This is the main script for performing reconstruction of ground contamination.

Required packages: mvtnorm and RColorBrewer

By default the script generates 10x10 activity matrix of 100 kBq/m^2 with one value of 350.
![image](https://user-images.githubusercontent.com/30175279/129942429-aba1484d-1de2-478a-b851-b66919175b0d.png)


Then, detector response is evaluated.

Using the evaluated response, estimated activity values are reconstructed.

![image](https://user-images.githubusercontent.com/30175279/129942499-5178280f-e6a2-446d-8085-c80108b5b470.png)


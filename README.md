Traffic sign detection in traffic images: the main purpose of this project is to detect and classify the traffic road signs in different traffic images. This project aims to create a transparent software system whose input is an image captured by a digital camera mounted on a vehicle. 
It is only available for European traffic signs.
For each image, the system returns a list of bounding-box coordinates and the specific class for every traffic sign detected in the image. Based on color and shape, the system is capable of detecting 7 main classes: STOP signs, Priority signs, Yield signs, Mandatory signs, Guidance/Information signs, Warning signs and Prohibitory signs. 
The purpose of this system is not to be so fast that it could work in real-time, but to explore many image processing techniques and computer vision algorithms for color-based segmentation, clustering and shape recognition. 
All the methods used in this project can be seen, analyzed, understood and controlled by the user, so that the AI system created can be transparent. 
The system has been tested using a dataset which consists of 379 traffic images manually labeled for this specific classes: they belong to public traffic sign recognition datasets such as ,,Road Signs Dataset” and ,,DFG Traffic Sign Data Set“ and public images on the internet. 
The overall precision is 88.45%. The overall recall is 89.12%. 
The entire code is written in Matlab.

Contents of posted code sources:
- TSD_MAIN - the main code whose input is a traffic image and returns the specific bounding-boxes and classed detected - in order not to be so long, it uses a function 'TSD_function' which computes all the algorithms and methods used - acts like an interface with the user
- TSD_function - a function which includes the proposed method and algorithms for detecting traffic signs in an image - receives a rgb_image and returns a list of bboxes
- myDBSCAN - personal implementation of DB-scan clustering method used for extracting individual objects in a binary image for later analysis and shape recognition
- myGaussianfiler - personal implementation of gaussian filter for smoothing
- myFRS - Fast Radial Symmetry algorithm used for circle recognition
- myHarrisDet - Harris Corner Detector algorithm used in checkCorners function
- checkCorners - a function whose input is a binary image containing one shape and checks specific zones if they are corners or not - the corner decision decides the shape
- myEOmeasures - geometrical shape measures for decision circle/octagon
- check...Area - Yellow/White/Blue - checks the inside color of the analyzed object

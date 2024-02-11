# Spicule-Oscillations
MATLAB script that finds physical characteristics (e.g. oscillation frequency, size) of solar spicules. It can also work for spicule-like structures such as tall grass when fine-tuned. 



## Using the spicule detection and analysis Matlab script

  
Note: This script was written in Matlab 2019b. It may not function properly when run in other versions.

1.  Setting up the script.

    1.  This script requires the Matlab Image Processing Toolbox which could be downloaded from “Apps" \> “Get more apps" \> search for “Image Processing Toolbox" \> add it to Matlab.

    2.  Download the folder containing the main Matlab script and other functions necessary for it to function.

    3.  Obtain a video of spicules or spicule-like structures to be detected and analysed, then copy it to the extracted folder.

2.  Setting up the main script settings.

    1.  Open the main script named “Main.m" and change the value of the “videoName" variable to the name of the video to be analysed. For example, “spicules.mov" where “spicules" is the name of the video and “.mov" indicates its format.

    2.  If you wish to check the input video, move to the "Play original frames" section and input values for the “vid.CurrentTime" and “playSpeed" variables corresponding to the start time and play speed of the video, respectively, and running with CTRL and ENTER. Similarly, you could check individual frames from the “Extract single frame" section by inputting a time value for the “vid.CurrentTime" variable and pressing CTRL and ENTER.

    3.  Set the “numFrames" variable’s value to 1 for the initial setup and adjustment of parameters.

    4.  Enter a value for the “realFrameWidth" variable which represents the real width in meters of the video to be analysed.

    5.  Enter a value for the “framePerSecond" variable which represents the time in seconds between each frame. If the video was captured in real-time, then enter the value as “vid.FrameRate".

    6.  Set the “newCrop" value to 1 to select a new region of interest when running the script. A value of 0 utilises the previous selected region of interest (if it exists).

3.  Following this step, it is advisable to regularly run the script using CTRL and ENTER to check for the effects of the choice of parameters.

    1.  Select the region of interest by clicking on the frame and dragging to create a box. This box could be readjusted and moved. Once the box encompasses the desired region, double-click on it to save its coordinates and proceed.

    2.  Create a new crop to exclude any detection artefacts.

    3.  Select the allowed dimension through which an object’s centroid could be tracked. This region will be used as a boundary to compare objects between frames and whether or not to select them as potential matches.

4.  Method of enhancement. Change the “enhanceMethod" variable to either 0, 1, 2, or 3, representing no image enhancement, “imadjust", “histeq", and “adapthisteq", respectively.

5.  Gray thresholding.

    1.  Change the “threshMethod" variable to 0, 1, or 2, corresponding to no thresholding, manual thresholding, and automatic thresholding.

    2.  If 1 was entered (corresponding to manual thresholding) then enter a value for “grayThreshold", otherwise ignore this step. A higher value raises the gray intensity threshold, allowing only dark objects, and a lower value allows both bright and relatively darker objects.

6.  Method of object detection. Enter a value for the “detectionMethod" variable where 1 performs edge detection while 2 performs curvilinear detection.

    1.  If 1 (corresponding to edge detection) was chosen

        1.  Enter a value for the “edgeMethod" variable between quotation marks. Choose one of the following: “sobel", “prewitt", “roberts", “log", “zerocross", “canny", or “approxcanny"

        2.  Enter either 0 for automatic parameters or 1 for manual parameters.

        3.  if 1 (corresponding to manual parameters) was selected, set the “edgeThresh" variable to a value between 0 and 1. Lower values increase the detection sensitivity.

        4.  If either the “sobel" or “prewitt" methods were chosen, set the “sobPrewDirect" to one of “vertical", “horizontal", or “both", corresponding to the direction of edges to be detected. Note: “vertical" covers 45 to 135 degree oriented objects, while “horizontal" includes orientations between 45 to -45 degrees.

        5.  If either the “log" or “canny" methods were chosen, adjust the “logCanSigma" variable. Higher values increase the grouping of edges and vice versa.

    2.  If 2 (corresponding to curvilinear detection) was chosen

        1.  Set a value of the “StepFilterGrammes" variable. Higher values join structures longitudinally while lower values cause greater separation in terms of length.

        2.  Set the “lowLimit" variable to a value between -6 and 0, corresponding to 180 to 90 degrees. Similarly, set the “highLimit" to a value between 0 and 6, corresponding to 90 to 0 degrees. These two variable set the range of orientations to be detected. The “StepFilterStep" variable controls the steps between the low and high boundaries.

        3.  Choose the step filter type by setting a value for the “StepFilterType" variable of either 0 for a polynomial filter or 1 for a step function.

        4.  Set the “BinThresh" variable to control the detection intensities allowed to remain. A lower value allows less detections of lower confidence (indicated by their intensity).

        5.  Set the “smallThresh" to a pixel value which eliminates detections with an area below that limit.

7.  Tracking the spicules.

    1.  Set the “newCrop2" variable to 0 to use the previously saved crop for artefact elimination. Set to 1 to use a new crop.

    2.  Similarly, set the “newTrackReg" to 0 to use the previous tracking region, and to 1 to use a new tracking region.

    3.  Set the “orientationLimit" variable to a value which represents the range of degrees to either side of an object where it could be considered the same object in the following frame.

    4.  Set the “matchIterations" variable to control the number of iterations to be completed while searching for a matching object between frames.

    5.  Set “showIndices" to 1 to display the tracking of the spicules. Setting to 0 does not plot the tracking process.

8.  Controlling which statistics to plot.

    1.  Set the “plotStats" variable to 1 to display the statistics of an individual chosen object.

    2.  Set the “plotCentOsc", “plotTopOsc", and “plotBottOsc" variables to 1 to display the displacement of the chosen object’s centroid, top, and bottom, respectively. Set to 0 to not display the displacement.

    3.  Set the “plotCentVelOsc", “plotTopVelOsc", and “plotBottVelOsc" to 1 to display the velocity of the chosen object’s centroid, top, and bottom, respectively. Set to 0 to not display the velocity.

    4.  Set the “plotDimensions" variable to 1 to display the length and width of the chosen object. Set to 0 to not display the dimensions.

    5.  Set the “plotOrientation" variable to 1 to display the orientation of the chosen object. Set to 0 to not display the orientation.

    6.  Set the “plotHistograms" variable to 1 to display all the detected objects’ statistics, or to 0 to not plot the histograms.

        1.  If the histograms are to be plotted, the “binNumber" variable controls the number of bins to be displayed in the histograms.

        2.  If the histograms are to be plotted, set the “numMatches" variable to a value between 1 and 4, indicating the top n distribution matches for the data displayed.

9.  After adjusting all the parameters, set the “numFrames" variable to the number of frames to be processed, and run the script.

10. Complete any required cropping operations (these will vary depending on the chosen settings). Note: their purpose will be described through text outputs in the Matlab console.

11. A new message will appear in the console, prompting an input to select an object to be analysed by its index. Ignore this step if the chosen settings do not request individual plotting. <span id="individSpic" label="individSpic"></span>

12. If individual plotting was requested, frame-by-frame plots will be constructed of the object’s top, centroid, bottom, their displacement and velocity (and any oscillation), length, width, and orientation. <span id="plotIndSpic" label="plotIndSpic"></span>

13. Messages concerning the existence of oscillations in certain object properties will appear in the console and histograms detailing overall object statistics will be displayed. <span id="rerun" label="rerun"></span>

14. A new message will appear in the console, prompting an input to either end the script with 0, or restart the script to plot an individual object with 1. If 1 was chosen to rerun the individual spicule detection, parts <a href="#individSpic" data-reference-type="ref" data-reference="individSpic">[individSpic]</a>, <a href="#plotIndSpic" data-reference-type="ref" data-reference="plotIndSpic">[plotIndSpic]</a>, and <a href="#rerun" data-reference-type="ref" data-reference="rerun">[rerun]</a> will repeat.

15. If you wish to view the binary frames that were analysed, visit the “Play uint8 frames" section, input loop and play speed values for the “loopLimit" and “playSpeed" variables, respectively, then press CTRL and ENTER. Similarly, a single binary frame could be played by visiting the “Play nth uint8 frame" section, inputting the frame to be displayed for the “displayFrame" variable, and pressing CTRL and ENTER.

16. If you have chosen the curvilinear detection method, you could view the detected frames by visiting the “Play double type frames" section, input loop and play speed values for the “loopLimit" and “playSpeed" variables, respectively, then press CTRL and ENTER. Similarly, a single double type frame could be played by visiting the “Play nth double type frame" section, inputting the frame to be displayed for the “displayFrame" variable, and pressing CTRL and ENTER.

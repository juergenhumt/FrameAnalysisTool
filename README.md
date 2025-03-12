# FrameAnalysisTool
Program for semi automatich approximation of higly complex ship or aircraft frames.

Using biarc approximation for frame representation of ship lines

A mathematical description of the frames of a wide variety of hull forms is challenging, because of, sometimes, not only numerous changes of radii of curvature along the frame, but a reversed direction of curvature as well. But such a description is needed for CAD models of ship hulls. Fortunately the proper tools for the job have been developed quite some time ago, in the 1960ies of the last century, as far as I know. The technique is called biarc approximation which means that the frames are divided into sections which are then approximated by two arc segments. Using two arc segments allows to connect points along the frame prescribing also the tangent at the start and endpoint. This allows to generate a smooth contour along even the most tortuous of frames. FrameAnalysisTool (= FAT) uses this technique. To analyze your own frames you have to generate an input file with pairs of coordinates, either in spread sheet csv format or in SVG format. See FAT descripton file

# Masterclass 22.09: Using path collective variables to find reaction mechanisms in complex free energy landscapes

This lesson was given as part of the PLUMED masterclass series in 2022.  It includes:

* Two videos that describe the theory. 
* A series of exercises for you to complete.
* A python notebook that contains a full solution to the exercises.

The flow chart shown below indicates the order in which you should consult the resources.  You can click on the nodes to access the various resources.  Follow the thick black lines for the best results.  The resources that are connected by dashed lines are supplmentary resources that you may find useful when completing the exercise.

This lesson was the eigth masterclass in the 2022 series.

```mermaid
flowchart TB;
  A[ref1] -.-> B[Lecture I];
  A -.-> C[Slides];
  B ==> D[Instructions];
  C -.-> D;
  D ==> E[Lecture II];
  E ==> F[notebook];
  click A "ref1" "A previous tutorial that introduces the basics of PLUMED syntax";
  click B "video1" "A lecture that was given on June 6th 2022 as part of the plumed masterclass series that introduces you to the exercises in this lesson";
  click C "masterclass_slides_s.pdf" "The slides that were used in the lecture on the 6th June";
  click D "INSTRUCTIONS.md" "Instructions for the exercise that you are supposed to complete";
  click E "video2" "A lecture that was given on June 13th 2022 as part of the plumed masterclass series that goes through the solutions to the exercises in the lesson";
  click F "exercises/masterclass-22-09.ipynb" "A python notebook that shows you how you can analyse some of the simulation data you generated";
```

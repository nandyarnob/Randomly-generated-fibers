# Randomly-generated-fibers
The program is used to randomly generate fibers between a cell and boundary. The co-ordinates of each of the fiber and the intersecting points are stored in various files. 

Here 36000 fibers are generated between the area of a circular cell and a square boundary. The cell radius is 7.6 mm and each fiber length is 7 mm. The boundary square has each side 30 times bigger than the radius of the cell. 

Here I used uniform random distribution module from Numpy and generated a single point which will be the middle point of the fiber. I also used the same random module to generate a corresponding slope angle. Using the length of the fiber, slope, and the middle point, I got the co-ordinates of  2 end points. I had to make sure that the middle point must be inside the boundary square and outside the cell circle. Using a while loop I generated required number of fibers within the constrains.

For using in ABACUS, I also needed the co-ordinates of the intersection of the 2 fibers. So I used simple algebra to get the co-ordinate of the intersection points. I also required the intersection point of any fiber and the cell circle. This was slightly more difficult. However simple quadratic equation was needed to solved to get the desired co-ordinates. 

After getting the intersecting points, I trimmed the isolated end points which are not part of any node. This ensures the proper generation of network of fibers and also be suitable to be used in ABACUS simulation. 

Finally, these intersecting co-ordinates are saved in different files as per ABACUS requirement. 

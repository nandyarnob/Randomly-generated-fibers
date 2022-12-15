# Randomly-generated-fibers
The program is used to randomly generate fibers between a cell and its boundary. The coordinates of each fiber and the intersecting points are stored in various files.

Here, 36,000 fibers are generated between the area of a circular cell and a square boundary. The cell radius is 7.6 mm, and each fiber length is 7 mm. The boundary square has each side 30 times bigger than the radius of the cell.

Here I used the uniform random distribution module from Numpy and generated a single point, which will be the middle point of the fiber. I also used the same random module to generate a corresponding slope angle. Using the length of the fiber, the slope, and the middle point, I got the coordinates of 2 end points. I had to make sure that the middle point was inside the boundary square and outside the cell circle. Using a while loop, I generated the required number of fibers within the constraints.

For use in ABACUS, I also needed the coordinates of the intersection of the 2 fibers. So I used simple algebra to get the coordinates of the intersection points. I also required the intersection point of any fiber and the cell circle. This was slightly more difficult. To obtain the desired coordinates, however, a simple quadratic equation had to be solved.

After getting the intersecting points, I trimmed the isolated end points, which are not part of any node. This ensures the proper generation of a network of fibers and is also suitable for use in ABACUS simulation.

Finally, these intersecting coordinates are saved in different files as per the ABACUS requirement.

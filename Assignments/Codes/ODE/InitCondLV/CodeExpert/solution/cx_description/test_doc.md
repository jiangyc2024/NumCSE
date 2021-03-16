There are three tests for `LV.hpp`, one for each variable printed in the following:
```
  std::pair<Vector2d, Matrix2d> PaW = PhiAndW(2.8, 1.5, 2);
  
  std::cout << "Test of PhiAndW():\nPhi = "
            << PaW.first.transpose()
            << "\nW = \n"
            << PaW.second << "\n\n";
  
  Vector2d y = findInitCond();
  std::cout << "Test of findInitCond():\ny = "
            << y.transpose() << "\n";
```
The correct output is
```
Test of PhiAndW():
Phi = 0.194835  2.25777
W = 
-0.158907 0.0637957
-0.121175 -0.610455

Test of findInitCond():
y = 3.10988 2.08098
```
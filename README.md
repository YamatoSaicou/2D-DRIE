# 2D-DRIE
This project is my graduation design in southeast university. It uses the cellular automata algorithm to convert the reactive ion etching system into the cellular automata system and constructs the etching physical model and follows the etching process. 
The project is based on C++ due to the professor's require. The result of the simulation is given out on the text.txt. And it is two-dimensional matrix, you can use matlab to draw it.
## The matlab draw script
clc
A=load('result.txt');
for y=50:100
    for x=1:50 
        switch A(x,y)
            case 2
                rectangle('Position',[y,-x,1,1],'facecolor','g')
                rectangle('Position',[100-y,-x,1,1],'facecolor','g')
            case 3
                rectangle('Position',[y,-x,1,1],'facecolor','w')
                rectangle('Position',[100-y,-x,1,1],'facecolor','w')
            case 4
                rectangle('Position',[y,-x,1,1],'facecolor','b')
                rectangle('Position',[100-y,-x,1,1],'facecolor','b')
            case 5
                rectangle('Position',[y,-x,1,1],'facecolor','b')
                rectangle('Position',[100-y,-x,1,1],'facecolor','b')
        end
    end
end
                
        
    

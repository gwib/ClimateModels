%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % INTRODUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % This exercise script is created as a basic introduction to people 
% % completely, or almost, new to using MATLAB. It therefore takes of 
% % from absolute scratch with increasing complexity into using some of the
% % calculations, functions and mindsets used in the course to which it is
% % created namely the joint DTU/Univ. of CPH course:
% % (KU): Climate?Models?and?Observations
% %?(DTU):?Climate?Models,?Observations?of?the?Past?and?the?Present,?
% % and?Projected?Climate?Change?including?Sea?Level?Rise?2018
% % 
% % The script is created on a windows based machine running the 2017a
% % MATLAB version. Using a Mac version and/or using another MATLAB release
% % (usually updated twice a year entitled with the year and either "a" or
% % "b"), might require minor changes to the script.
% % 
% % In the following an exercise begins with an !!
% % A new topic begins is marked with: %%%%%%%%%%%%%%%%%
% % 
% % Good luck and best regards, Morten Andreas Dahl Larsen, DTU
% % Contact me if you have questions or corrections: madla@dtu.dk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New script %%%%%%%%%%%%%%%%%%%%

% Start new script (Ctrl+N) or open existing script (Ctrl+O) from your machine 
% in the HOME or EDITOR tabs. Save scripts in the EDITOR tab (Ctrl+S).
% If you close MATLAB with a script (or more) in the Editor window it will load with this same script(s) when MATLAB is started
% again.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to run a script?  %%%%%%%%%%

% 'Selected possibilites often used by the author'
% 
% 1) You run the entire script by pressing F5 or pressing "Run" in the
% "EDITOR" tab above.
% 2) You can a section of the script by highlighting it and pressing F9.
% 3) You can divide your script into sections and then pressing Ctrl+ENTER
% to run the section in which the cursor is placed.
%       A section break is defined by "%%" (without spacing in between).
%       Sections are useful if you're developing code and only want to run 
%       the part in which you're working. Also if previous parts of the script
%       are 'heavy', taking a longer time to run - in the end it's all
%       about creating effective code, not least in terms of computation time
%       spent.
%       !! Try to implement sections when you continue with the exercise
%       !! below

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear workspace %%%%%%%%%%%%%%%

% You can start the running of a script with this statement since previously 
% used variables (as visible in the workspace window) can still be loaded and 
% affect the current run (if names overlap)
clear all

% !! When you do a script/exercise further down, try to do it with/without
% !! "clear all" to see the difference in the workspace window

% You can also clear single variables. We here define three variables:

A = 1;
AA = 1.5;
B = 2;

% And clear one of them:

clear A

% Or clear the ones starting with "A" after redefining A=1:

A = 1; % see that it reappears in the workspace

% And then:

clear A* % now only the "B" variable is left

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order and directories %%%%%%%%%%

% It is important to keep track of your directories and have a good and 'clean' 
% subdirectory system for (all of) your work. Therefore;
% !! Save this script in the a designated folder on your computer. Not on
% the desktop.
% Thereafter, it is equally important to define your working directory in
% MATLAB - this means that the script and files (including plots) created
% in this relation is saved in the same, and correct, folder.
% !! 1) Define the correct working directory in the directory window above in
% !! the MATLAB GUI (graphical user interface).
% !! Also try to do it using the "cd" function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Outcommenting" %%%%%%%%%%%%%%%%

% Ctrl+R outcomments a statement making it inactive (if you're 1) not
% using it, 2) still developing it, 3) if you need to write text to keep
% track of your script (as this line))
% Ctrl+T does the opposite
% You can also just manually write "%" in the beginning of a line

% !! Try to use these commands!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple equations %%%%%%%%%%%%%%%

% Simple math equations work just like you would write them on a piece of
% paper.
% Examples are given below

2+2;
9*10;
17-8;
70/7;

% However, the above equations do not define new variables, and therefore 
% only the latest is saved in the workspace window - entitled "ans". 
% !! Check for yourself.
% !! Try to define new variables based on these (simple) calculations above
% !! by writing the variable name followed by "=" in front of the equation,
% e.g: 

Course_years = 2018-2014+1; % instead of just 
2018-2014+1;

% !! Instead of writing the above years in numbers these variables could
% !! also have been defined prior to the calculation. This is often used and
% !! creates a flexibility to change the outcome of the calculations by
% !! criteria you define in the beginning of the script:

Current_year = 2018;
Course_start_year = 2014;
Course_years2 = Current_year - Course_start_year + 1;

% OBS - variables 1) cannot start with a number, 2) are case-sensitive and
% 3) must not include "-", use instead "_" (to avoid mixups with the
% subtraction sign)

% !! Try to create you own simple math equations, also defining new
% variables - be creative.
a = 6;
b = 37;

f = a*b;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices %%%%%%%%%%%%%%%%%%%%%%%

% MATLAB is a powerful tool to work with matrices whether in 1,2,3 or 4 dimensions 
% (or even further) - hence the name 
% A 1D matric (or vector) is defined by hard brackets, e.g.:

ONE_D_example_row = [0	0	0	0	0	0	60	322	607	829	847	923	856	734	791	590	319	124	27	0	0	0	0	0];
% (These data are hourly global radiation (W/m2) for a summer day in sourthern 
% Italy) 

% Using semicolons changes the row so that in the example below each number 
% is written to a new row.:

ONE_D_example_column = [0;0;0;0;0;0;60;322;607;829;847;923;856;734;791;590;319;124;27;0;0;0;0;0];

% !! See the difference in the workspace.

% Using "'" (end of matrix) transposes the matrix (to give the same alignment as 
% the above):

ONE_D_example_row_transposed = [0	0	0	0	0	0	60	322	607	829	847	923	856	734	791	590	319	124	27	0	0	0	0	0]';

% As for single  variables or unique number, matrices can also part of
% equations. 

% Examples:

% You want to multiply your time series by a factor of 1.25:
ONE_D_example_multiply_single_factor = ONE_D_example_column * 1.25; % check that it works and try your own example

% You want to multiply by a varying time factor throughout the 24h.
% First create the time series to multiply with:
Factor = rand(24,1)+0.5; % creates a 1*24 matrix of random numbers between 0.5 and 1.5
% Plot the multiplication matrix 
figure
plot(Factor)
% Perform the multiplication:
ONE_D_example_multiply_varying_factor = ONE_D_example_column.*Factor;
% Be aware here that we need to add the "." before "*" to tell MATLAB that we
% multiply one-by-one

% Now lets plot the original time series as well as the two modified ones:
figure % build figure
hold on % tell MATLAB that the forthcoming plot calls should be on the same plot
plot(ONE_D_example_column)
plot(ONE_D_example_multiply_single_factor)
plot(ONE_D_example_multiply_varying_factor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Documentation, help % File Exchange %%%%%%%%

% Matlab functions and specifics are widely documented, e.g. 
% http://se.mathworks.com/help/matlab/

% You can also right-click on a function and press "help on selecetion" or press F1
% after highlighting a specific function

% Another common help source for more advanced questions include Stackoverflow:
% http://stackoverflow.com/questions/tagged/matlab

% It is also possible to create your own, more universal, functions which are ".m" 
% script based just like the present script (.m is the MATLAB script format but
% can be opened in any text editor). The library for user created functions of
% this sort (File Exchange) is located here: 
% http://www.mathworks.com/matlabcentral/fileexchange/, and can be a useful resource
% to do things not implemented in MATLAB in the downloadable version.
% Notice that the functions are 'graded' by users which could be a
% guideline for whether its usability and correctness is appropriate.
% When downloading external functions they need to be placed in the appropriate 
% folder and be linked to. E.g.:
% 'C:\Program
% Files\MATLAB\R2015b\toolbox\NameOfFunction\...'
% Hereafter, go to the "HOME" tab, "Set Path", add that folder and press "Save"

% Assignments for the above:
% 1) !! go to "help" for the "plot" function used previously and change e.g. the curve
% colors, the curve linewidth and the marker/line style
% 2) !! go to the "help" section for the function "repmat" and create a 2D
% matrix of your own choice
% 3) !! download a function of your own choice i File Exchange and use it
% in an example of your own choice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppress output and time consumption %%%%%%%%%

% All of the above calculations are, as default, written to the command window 
% below making you able to follow the progress.
% However, it also makes the calculations considerably slower when
% involving more extensive data sets. You can therefore suppress output to
% the command window by ending expression with i ";".

% Further, you can time your calculations by using the commands "tic" and "toc"
% before/after what you want to time

% Here's an example highlightning both of the above:

% Make a data set (10x10) of random numbers:
tic
rand(10);
Without_suppressing_10 = toc;

tic
rand(10);
With_suppressing_10 = toc;

% Now, let's make a larger data set of 100x100 random numbers:
tic
rand(1e2);
Without_suppressing_100 = toc;
tic
rand(1e2);
With_suppressing_100 = toc;

% Let's try with 500x500:
tic
rand(500);
Without_suppressing_500 = toc;
tic
rand(500); 
With_suppressing_500 = toc;

% Let's try with 2000x2000:
tic
rand(2e3);
Without_suppressing_2000 = toc;
tic
rand(2e3); 
With_suppressing_2000 = toc;

% That's a factor of:
Speed_factor_10 = Without_suppressing_10/With_suppressing_10;
Speed_factor_100 = Without_suppressing_100/With_suppressing_100;
Speed_factor_500 = Without_suppressing_500/With_suppressing_500;
Speed_factor_2000 = Without_suppressing_2000/With_suppressing_2000;

% Let's plot the results:
figure
plot([10,100,500,2000],[Speed_factor_10,Speed_factor_100,Speed_factor_500,Speed_factor_2000],'*--k')

% You can see where this is going, and why I stopped at 2000x2000
% With 3D data the differences would have been even larger

% Notice, you can also write calculations in the command window as a simple 
% calculator, although these % are not saved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual appearance of code %%%%%%%%%%%%%%%%%%%%

% You can have more expressions on one line by seperating with a ";":
ABC = 7*7; DEF = 9*9; 

% You can split a long command into several lines by using "..."
figure
plot([10,100,1000],[Speed_factor_10,Speed_factor_100,Speed_factor_2000],...
    'color',[1 0.2 0.2],...
    'linewidth',1.5,...
    'marker','o',...
    'markersize',12,...
    'MarkerEdgeColor',[0.2 0.2 1],...
    'MarkerFaceColor','g')

% Spacing in expressions does not matter, so choose the style you like the most:

HIJ1 = ABC + DEF;
HIJ2=ABC+DEF;

% Work with logical variable names, i.e.
%     1) Not too long (lack of overview)
%     2) Not too short (difficulty in knowing what they are)
%     3) Think about starting letter, so they group according to your preference 
%     in the workspace window
% 
% And; I likely did not succed in living up to any of the above rules in this 
% exercise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This exercise already used the 'plot' function.
% Other common plotting functions include:

% Data to plot - import data from an external Excel sheet:
Plot_data = xlsread('Matlab_intro_data.xlsx','Ea'); % only works if the Excel sheet is located in the
% same folder as this script. Otherwise, use the entire directory 'C:\...\...'
% The data are simulated daily actual evapotranspiration (mm) over Southern 
% Italy in 2012.
% And the corresponding x-axis data is, when necessary, just a series from 1 to 
% 365 (days): 
Plot_time = [1:365];

% !! Try to change the specifications for each plot

% 1) plot
figure
plot(Plot_time,Plot_data(:,2));

% 2) scatter
figure
scatter(Plot_time,Plot_data(:,2));

% 3) bar
figure
bar(Plot_data(:,2));

% 4) histogram
edges = [0:0.25:5];
figure
histogram(Plot_data(:,2),edges);

% And for distributed (2D) plotting I like imagesc. First, load the data;
% in this case average yearly precipitation (mm) from the E-OBS data set in 
% the period 1950-2014:

Plot_data_dist = xlsread('Matlab_intro_data.xlsx','Precip');
Plot_data_dist = flipud(Plot_data_dist); % flip up/down

% 5) imagesc
figure
hold on
imagesc(Plot_data_dist)
caxis([100 3000]) % set the coloraxis min and max
colormap(jet) % colormap

% If you want to suppress the figure output (while e.g. still writing to a file)
% use figure('visible','off')

% Save a figure to file by using e.g. print 

figure('visible','off')
hold on
imagesc(Plot_data_dist)
caxis([100 3000]) % set the coloraxis min and max
colormap(jet) % colormap
print('-dtiff','-r300','PRECIP_plot_test'); % also stating resolution and format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra plotting functions %%%%%%%%%%%%%

% Other nice functions here exemplified to the "plot" function:

%%

figure('visible','off')
plot(Plot_time,Plot_data,'b','linewidth',1.1); % same as above
set(gca,'XTick',[0:365.25/12:366]); % define x-tick placement, here roughly every month
set(gca,'XTickLabel',{'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'}); % x-tick labels
set(gca,'YTick',[0:1:5]); % define y-tick placement
ylim([0 5]); % y-axis min and max
xlim([0 366]); % x-axis min and max
set(gca,'Fontsize',10); % font size
xlabel('Date','Fontsize',11,'Fontweight','bold');
ylabel('Ea (mm/day)','Fontsize',11,'Fontweight','bold');
title('Actual daily evapotranspiration, Southern Italy, 2012','Fontsize',12)
set(gca,'LineWidth',1)
text(177,1,'Rainy period?','color','r','fontsize',12) % text example
box on % figure edges
grid off % internal grid
set(gcf, 'Color', 'w'); % graph color
print('-dtiff','-r300','Ea no legend'); % save without legend

h_legend = legend('Ea (mm/day)','Location','NorthWest'); % Legend text and location
set(h_legend,'Fontsize',11,'Fontweight','bold'); % Legend characteristics

print('-dtiff','-r300','Ea with legend'); % save WITH legend
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For loops %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For loops are needed to perform, typically larger and more complex, operations 
% a number of times dependent on e.g. % a variable that you then modify for each 
% loop

% In this example we will take the mean for every 3x3 cell in our EOBS precipitation
% data set - i.e. an upscaling.


Precip_EOBS = Plot_data_dist; % Let's rename the data no so we know that we are looking at precipitation

% These are x and y positions for the new upscaled file and they need to be 
% defined as zero before the loop
x_upscale = 0; 
y_upscale = 0;


% It is a good idea clearing the resulting variable before going into a loop - 
% It might appear from previous development work 
clear Precip_EOBS_upscaled

% We start the by stating that we want to go through x values from 1 to the lenght of 
% the x-direction minus 2 in steps of 3
for x = 1:3:numel(Precip_EOBS(:,1))-2;
%     We then want x_upscale to increase by 1 each time the loop reaches
%     this position because here we know that the original x value is
%     increased
    x_upscale = x_upscale + 1;
%     The y_upscale need to be reset for every x-direction loop because we
%     start by looping over the x-direction and then run all the y's.
    y_upscale = 0;
%     We now loop in the y-direction for each x-value, still jumping 3
%     steps until reaching the lenght of the y-dimension minus 2 
    for y = 1:3:numel(Precip_EOBS(1,:))-2;
%         Now inside the y loop we increment by 1 for each loop
            y_upscale = y_upscale + 1;
%             We here calculate the new upscaled precipitation placing results in 
%             the positions of the x_upscale/y_upscale values which have one 'spot'
%             for every 3 original x/y spot
%             The mean is calculated from the x/y value in question
%             (depending on how far we are in the loop) and then 2 spots
%             ahead. We use "nanmean" because the data have "nan" values
%             meaning that there are no data here (ocean = no measurement stsation)
Precip_EOBS_upscaled(x_upscale,y_upscale) = nanmean(nanmean(Precip_EOBS(x:x+2,y:y+2)));
%     For each "for" there is an "end" closing the loop
    end
end


Precip_EOBS_upscaled = flipud(Precip_EOBS_upscaled); % flip the data so they are plotted correctly
figure % plot the data
imagesc(Precip_EOBS_upscaled)

% Compare with
Precip_EOBS = flipud(Precip_EOBS);
figure 
imagesc(Precip_EOBS)

% !! Understand what goes on above and try to change what it does. E.g.
% !! averaging for each 4th cell.


// This example macro demonstrates how to use the getProfile() function

  // Set alt key down for vertical profile of a rectangular selection
  //setKeyDown("alt");

  // Get profile and display values in "Results" window
  //setBatchMode(true);

  run("Clear Results");
  width = getWidth;
  height = getHeight;
  depth = nSlices;
  name = getTitle;
  lineCount = 0;
  actionToDo = "Continue";
  
  while(actionToDo == "Continue"){
      lineCount++;
	  waitForUser("Fluorescence profile macro","Draw a line to get profile");
	  for (k=1; k<=depth; k++) {
		  	profile = getProfile();
		  	setSlice(k);
			for (i=0; i<profile.length; i++) {
			  setResult("Profile_Slice_" + k, i, profile[i]);
			}
			updateResults;
	  }
	
	  path = getDirectory("C://Users//JosephVermeil//Desktop//Image Analysis//FluoProfile//") + name + "_Line_" + lineCount + ".txt";
	  saveAs("Results", path);
	
	  Dialog.create("Fluorescence profile macro");
	  items = newArray("Continue", "Stop");
	  Dialog.addRadioButtonGroup("Action", items, 2, 1, "Continue");
	  Dialog.show;
	  actionToDo = Dialog.getRadioButton;

  }
  
  //setBatchMode(false);
	
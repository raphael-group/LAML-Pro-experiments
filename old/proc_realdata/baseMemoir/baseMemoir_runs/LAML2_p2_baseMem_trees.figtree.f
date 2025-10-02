#NEXUS
begin taxa;
	dimensions ntax=30;
	taxlabels
	pos2_cell1
	pos2_cell10
	pos2_cell11
	pos2_cell12
	pos2_cell13[&!color=#ff0000]
	pos2_cell14
	pos2_cell15[&!color=#ff0000]
	pos2_cell16
	pos2_cell17[&!color=#ff0000]
	pos2_cell18[&!color=#ff0000]
	pos2_cell20
	pos2_cell23[&!color=#ff0000]
	pos2_cell24
	pos2_cell25[&!color=#ff0000]
	pos2_cell26[&!color=#ff0000]
	pos2_cell27[&!color=#ff0000]
	pos2_cell28[&!color=#ff0000]
	pos2_cell29
	pos2_cell3
	pos2_cell30
	pos2_cell31
	pos2_cell32
	pos2_cell33
	pos2_cell34
	pos2_cell35
	pos2_cell4
	pos2_cell5
	pos2_cell6
	pos2_cell8[&!color=#ff0000]
	pos2_cell9
;
end;

begin trees;
	tree tree_1 = [&R] ((((((((pos2_cell12:0.063482,pos2_cell3:0.063482):0.441719,(pos2_cell13:0.006008,pos2_cell6:0.006008):0.499193):0.421133,((pos2_cell25:0.669226,pos2_cell16:0.669226):0.232623,((pos2_cell24:0.485662,pos2_cell30:0.485662):0.412538,(((pos2_cell31:0.01006,pos2_cell4:0.01006):0.472564,pos2_cell1:0.482624):0.274856,pos2_cell23:0.75748):0.14072):0.003649):0.024485):5.95E-4,(((pos2_cell17:0.143886,((pos2_cell26:4.7E-4,pos2_cell28:4.7E-4):4.7E-4,pos2_cell27:9.39E-4):0.142947):0.425966,pos2_cell8:0.569852):0.325877,pos2_cell29:0.895729):0.031199):0.011405,pos2_cell34:0.938333):0.021649,((pos2_cell20:0.662931,pos2_cell11:0.662931):0.225834,(((pos2_cell18:0.427357,pos2_cell32:0.427357):0.411225,(pos2_cell33:0.194535,pos2_cell35:0.194535):0.644047):0.04829,(pos2_cell5:0.202902,pos2_cell10:0.202902):0.683971):0.001893):0.071216):0.017515,((pos2_cell14:0.232887,pos2_cell9:0.232887):0.19554,pos2_cell15:0.428426):0.549071):0.022503);
end;

begin figtree;
	set appearance.backgroundColorAttribute="Default";
	set appearance.backgroundColour=#ffffff;
	set appearance.branchColorAttribute="User selection";
	set appearance.branchColorGradient=false;
	set appearance.branchLineWidth=1.0;
	set appearance.branchMinLineWidth=0.0;
	set appearance.branchWidthAttribute="Fixed";
	set appearance.foregroundColour=#000000;
	set appearance.hilightingGradient=false;
	set appearance.selectionColour=#2d3680;
	set branchLabels.colorAttribute="User selection";
	set branchLabels.displayAttribute="Branch times";
	set branchLabels.fontName="sansserif";
	set branchLabels.fontSize=8;
	set branchLabels.fontStyle=0;
	set branchLabels.isShown=false;
	set branchLabels.significantDigits=4;
	set layout.expansion=0;
	set layout.layoutType="RADIAL";
	set layout.zoom=0;
	set legend.attribute=null;
	set legend.fontSize=10.0;
	set legend.isShown=false;
	set legend.significantDigits=4;
	set nodeBars.barWidth=4.0;
	set nodeBars.displayAttribute=null;
	set nodeBars.isShown=false;
	set nodeLabels.colorAttribute="User selection";
	set nodeLabels.displayAttribute="Node ages";
	set nodeLabels.fontName="sansserif";
	set nodeLabels.fontSize=8;
	set nodeLabels.fontStyle=0;
	set nodeLabels.isShown=false;
	set nodeLabels.significantDigits=4;
	set nodeShapeExternal.colourAttribute="User selection";
	set nodeShapeExternal.isShown=false;
	set nodeShapeExternal.minSize=10.0;
	set nodeShapeExternal.scaleType=Width;
	set nodeShapeExternal.shapeType=Circle;
	set nodeShapeExternal.size=4.0;
	set nodeShapeExternal.sizeAttribute="Fixed";
	set nodeShapeInternal.colourAttribute="User selection";
	set nodeShapeInternal.isShown=false;
	set nodeShapeInternal.minSize=10.0;
	set nodeShapeInternal.scaleType=Width;
	set nodeShapeInternal.shapeType=Circle;
	set nodeShapeInternal.size=4.0;
	set nodeShapeInternal.sizeAttribute="Fixed";
	set polarLayout.alignTipLabels=false;
	set polarLayout.angularRange=0;
	set polarLayout.rootAngle=0;
	set polarLayout.rootLength=100;
	set polarLayout.showRoot=true;
	set radialLayout.spread=0.0;
	set rectilinearLayout.alignTipLabels=false;
	set rectilinearLayout.curvature=0;
	set rectilinearLayout.rootLength=100;
	set scale.offsetAge=0.0;
	set scale.rootAge=1.0;
	set scale.scaleFactor=1.0;
	set scale.scaleRoot=false;
	set scaleAxis.automaticScale=true;
	set scaleAxis.fontSize=8.0;
	set scaleAxis.isShown=false;
	set scaleAxis.lineWidth=1.0;
	set scaleAxis.majorTicks=1.0;
	set scaleAxis.minorTicks=0.5;
	set scaleAxis.origin=0.0;
	set scaleAxis.reverseAxis=false;
	set scaleAxis.showGrid=true;
	set scaleBar.automaticScale=true;
	set scaleBar.fontSize=10.0;
	set scaleBar.isShown=true;
	set scaleBar.lineWidth=1.0;
	set scaleBar.scaleRange=0.0;
	set tipLabels.colorAttribute="User selection";
	set tipLabels.displayAttribute="Names";
	set tipLabels.fontName="sansserif";
	set tipLabels.fontSize=15;
	set tipLabels.fontStyle=0;
	set tipLabels.isShown=true;
	set tipLabels.significantDigits=4;
	set trees.order=false;
	set trees.orderType="increasing";
	set trees.rooting=false;
	set trees.rootingType="User Selection";
	set trees.transform=false;
	set trees.transformType="cladogram";
end;


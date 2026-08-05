(* ::Package:: *)

(* ::Package:: *)
(**)


BeginPackage["mpmPlot`"];

(* ==================== \:63a5\:53e3\:4e0e Usage \:8bf4\:660e ==================== *)
loadMPMData::usage = 
  "loadMPMData[dirPath_String]\n" <>
  "\:3010\:8f93\:5165\:3011: dirPath (2D dump_p \:7c92\:5b50\:6570\:636e\:6587\:4ef6\:5939\:8def\:5f84)\:3002\n" <>
  "\:3010\:8f93\:51fa\:3011: 2D \:7c92\:5b50\:6570\:636e\:7ed3\:6784 Association \:5217\:8868\:3002";

loadGridData::usage = 
  "loadGridData[dirPath_String]\n" <>
  "\:3010\:8f93\:5165\:3011: dirPath (2D dump_g \:80cc\:666f\:7f51\:683c\:6570\:636e\:6587\:4ef6\:5939\:8def\:5f84)\:3002\n" <>
  "\:3010\:8f93\:51fa\:3011: 2D \:7f51\:683c\:6570\:636e\:7ed3\:6784 Association \:5217\:8868\:3002";

mPlotData::usage = 
  "mPlotData[pDataVar, opts___] \:6216 mPlotData[pDataVar, gDataVar, opts___]\n" <>
  "\:3010\:8f93\:5165\:3011: pDataVar (\:7c92\:5b50\:53d8\:91cf\:540d), gDataVar (\:53ef\:9009, \:7f51\:683c\:53d8\:91cf\:540d)\:3002\n" <>
  "\:3010\:9009\:9879\:3011: \"defaultVar\", \"pointSize\", ColorFunction, ImageSize, \"exportPath\"\:3002";

loadMPMData3D::usage = 
  "loadMPMData3D[dirPath_String]\n" <>
  "\:3010\:8f93\:5165\:3011: dirPath (3D dump_p \:7c92\:5b50\:6570\:636e\:6587\:4ef6\:5939\:8def\:5f84)\:3002\n" <>
  "\:3010\:8f93\:51fa\:3011: 3D \:7c92\:5b50\:6570\:636e\:7ed3\:6784 Association \:5217\:8868\:3002";

loadGridData3D::usage = 
  "loadGridData3D[dirPath_String]\n" <>
  "\:3010\:8f93\:5165\:3011: dirPath (3D dump_g \:80cc\:666f\:7f51\:683c\:6570\:636e\:6587\:4ef6\:5939\:8def\:5f84)\:3002\n" <>
  "\:3010\:8f93\:51fa\:3011: 3D \:7f51\:683c\:6570\:636e\:7ed3\:6784 Association \:5217\:8868\:3002";

mPlotData3D::usage = 
  "mPlotData3D[pDataVar, opts___] \:6216 mPlotData3D[pDataVar, gDataVar, opts___]\n" <>
  "\:3010\:8f93\:5165\:3011: pDataVar (\:7c92\:5b50\:53d8\:91cf\:540d), gDataVar (\:53ef\:9009, \:7f51\:683c\:53d8\:91cf\:540d)\:3002\n" <>
  "\:3010\:9009\:9879\:3011: \"defaultVar\", \"defaultView\", \"pointSize\", ColorFunction, ImageSize, \"exportPath\"\:3002";

Begin["`Private`"];

(* \:5f3a\:5236\:5168\:4f53\:4f20\:5740\:ff0c\:9632\:6b62UI\:5047\:6b7b *)
SetAttributes[{mPlotData, mPlotData3D}, HoldAll];

If[Kernels[] == {}, LaunchKernels[]];



(* ========================================================================= *)
(* 1. \:9ad8\:6027\:80fd\:989c\:8272\:5411\:91cf\:5316\:6620\:5c04\:4e0e\:56fe\:4f8b\:683c\:5f0f\:5316 (\:79d1\:5b66\:8ba1\:6570\:6cd5)                            *)
(* ========================================================================= *)
$colorTable = Developer`ToPackedArray[List @@@ Table[ColorData["Rainbow"][i/255.], {i, 0, 255}], Real];

getColorIndices = Compile[{{val, _Real}, {minV, _Real}, {maxV, _Real}},
  Module[{scaled, idx},
    scaled = (val - minV) / (maxV - minV + 1.0*^-12);
    idx = Round[scaled * 255.] + 1;
    If[idx < 1, 1, If[idx > 256, 256, idx]]
  ],
  RuntimeAttributes -> {Listable},
  RuntimeOptions -> "Speed"
];

getFastColors[vals_, minV_, maxV_] := $colorTable[[getColorIndices[vals, minV, maxV]]];

(* \:2705 \:5f3a\:5236\:56fa\:5b9a\:56fe\:4f8b\:683c\:5f0f\:ff1a4\:4f4d\:6709\:6548\:6570\:5b57\:ff0c\:81ea\:52a8\:8865\:96f6\:ff0c\:53733\:4f4d\:5c0f\:6570\:7684\:79d1\:5b66\:8ba1\:6570\:6cd5\:ff0c\:9632\:6b62\:8fb9\:6846\:6296\:52a8 *)
sciTicks[min_, max_] := Table[{v, ScientificForm[N[v], 4, NumberPadding -> {"", "0"}]}, {v, FindDivisions[{min, max}, 5]}];

(* ========================================================================= *)
(* 2. \:5065\:58ee\:7684\:6d41\:5f0f\:6570\:636e\:89e3\:6790                                                       *)
(* ========================================================================= *)
fastParseRobust[filepath_String, is3D_ : False] := Module[
  {rawText, sections, headerLine, dataStr, headers, rawTable, cleanMat, cols, assoc, coords},
  rawText = Import[filepath, "Text"];
  sections = StringSplit[rawText, "ITEM: ATOMS "];
  If[Length[sections] < 2, Return[<||>]];
  headerLine = First[StringSplit[sections[[2]], {"\r\n", "\n"}]];
  headers = StringSplit[headerLine];
  dataStr = StringDrop[sections[[2]], StringLength[headerLine] + 1];
  rawTable = ImportString[dataStr, "Table"];
  cleanMat = Developer`ToPackedArray[Select[rawTable, VectorQ[#, NumericQ] &], Real];
  If[Length[cleanMat] == 0, Return[<||>]];
  cols = Transpose[cleanMat];
  assoc = AssociationThread[headers, cols];
  coords = If[is3D, Transpose[{assoc["x"], assoc["y"], assoc["z"]}], Transpose[{assoc["x"], assoc["y"]}]];
  <|"Headers" -> headers, "Data" -> assoc, "Coords" -> coords|>
];

(* ========================================================================= *)
(* 3. \:6570\:636e\:52a0\:8f7d\:51fd\:6570 (\:7c92\:5b50 & \:7f51\:683c)                                               *)
(* ========================================================================= *)
loadMPMData[dataDir_String] := Module[{files, res},
  files = SortBy[FileNames["dump_p.*.LAMMPS", dataDir], ToExpression[StringExtract[FileNameTake[#], "." -> -2]] &];
  If[Length[files] == 0, MessageDialog["\:672a\:627e\:5230 2D \:7c92\:5b50 dump_p \:6587\:4ef6\:ff01"]; Return[$Failed]];
  PrintTemporary["\:26a1 \:5e76\:884c\:52a0\:8f7d 2D \:7c92\:5b50\:6570\:636e..."];
  res = ParallelMap[fastParseRobust[#, False] &, files, DistributedContexts -> Automatic];
  Print["\:2705 2D \:7c92\:5b50\:5bfc\:5165\:6210\:529f\:ff0c\:5171 ", Length[res = Select[res, # =!= <||> &]], " \:5e27\:ff01"]; res
];

loadMPMData3D[dataDir_String] := Module[{files, res},
  files = SortBy[FileNames["dump_p.*.LAMMPS", dataDir], ToExpression[StringExtract[FileNameTake[#], "." -> -2]] &];
  If[Length[files] == 0, MessageDialog["\:672a\:627e\:5230 3D \:7c92\:5b50 dump_p \:6587\:4ef6\:ff01"]; Return[$Failed]];
  PrintTemporary["\:26a1 \:5e76\:884c\:52a0\:8f7d 3D \:7c92\:5b50\:6570\:636e..."];
  res = ParallelMap[fastParseRobust[#, True] &, files, DistributedContexts -> Automatic];
  Print["\:2705 3D \:7c92\:5b50\:5bfc\:5165\:6210\:529f\:ff0c\:5171 ", Length[res = Select[res, # =!= <||> &]], " \:5e27\:ff01"]; res
];

loadGridData[dataDir_String] := Module[{files, res},
  files = SortBy[FileNames["dump_g.*.LAMMPS", dataDir], ToExpression[StringExtract[FileNameTake[#], "." -> -2]] &];
  If[Length[files] == 0, MessageDialog["\:672a\:627e\:5230 2D \:7f51\:683c dump_g \:6587\:4ef6\:ff01"]; Return[$Failed]];
  PrintTemporary["\:26a1 \:5e76\:884c\:52a0\:8f7d 2D \:7f51\:683c\:6570\:636e..."];
  res = ParallelMap[fastParseRobust[#, False] &, files, DistributedContexts -> Automatic];
  Print["\:2705 2D \:7f51\:683c\:5bfc\:5165\:6210\:529f\:ff0c\:5171 ", Length[res = Select[res, # =!= <||> &]], " \:5e27\:ff01"]; res
];

loadGridData3D[dataDir_String] := Module[{files, res},
  files = SortBy[FileNames["dump_g.*.LAMMPS", dataDir], ToExpression[StringExtract[FileNameTake[#], "." -> -2]] &];
  If[Length[files] == 0, MessageDialog["\:672a\:627e\:5230 3D \:7f51\:683c dump_g \:6587\:4ef6\:ff01"]; Return[$Failed]];
  PrintTemporary["\:26a1 \:5e76\:884c\:52a0\:8f7d 3D \:7f51\:683c\:6570\:636e..."];
  res = ParallelMap[fastParseRobust[#, True] &, files, DistributedContexts -> Automatic];
  Print["\:2705 3D \:7f51\:683c\:5bfc\:5165\:6210\:529f\:ff0c\:5171 ", Length[res = Select[res, # =!= <||> &]], " \:5e27\:ff01"]; res
];

(* ========================================================================= *)
(* 4. 2D \:7ed8\:56fe\:6a21\:5757 (\:4e25\:683c\:4e09\:884c\:5168\:82f1\:6587\:7d27\:51d1\:5e03\:5c40 + \:52a8\:6001\:523b\:5ea6 + \:5b9a\:5bbd\:6392\:7248)               *)
(* ========================================================================= *)
Options[mPlotData] = {"defaultVar" -> "seq", "pointSize" -> 0.005, "exportPath" -> "", ColorFunction -> "Rainbow", ImageSize -> Large};

mPlotData[pData_, gData_, opts : OptionsPattern[]] /; !MatchQ[Unevaluated[gData], _Rule | _RuleDelayed] := Module[
  {hasGrid, numFrames, allHeaders, vars, defVar, exportDir, imgSz},
  
  hasGrid = (Unevaluated[gData] =!= None);
  If[!MatchQ[Unevaluated[pData], _Symbol], MessageDialog["Please assign particle data to a variable first!"]; Return[$Failed]];
  If[hasGrid && !MatchQ[Unevaluated[gData], _Symbol], MessageDialog["Please assign grid data to a variable first!"]; Return[$Failed]];
  
  numFrames = If[hasGrid, Min[Length[pData], Length[gData]], Length[pData]];
  If[numFrames == 0, Return["Data is empty"]];
  
  allHeaders = DeleteDuplicates[Flatten[#["Headers"] & /@ pData]];
  If[hasGrid, allHeaders = DeleteDuplicates[Join[allHeaders, Flatten[#["Headers"] & /@ gData]]]];
  vars = Select[allHeaders, ! MemberQ[{"id", "type", "tag", "x", "y"}, #] &];
  
  defVar = OptionValue["defaultVar"]; If[! MemberQ[vars, defVar], defVar = First[vars]];
  exportDir = OptionValue["exportPath"]; imgSz = OptionValue[ImageSize];

  Manipulate[
    Module[{pCoords, pField, gCoords, gField, minV, maxV, plotElements = {}},
      minV = Infinity; maxV = -Infinity;
      
      If[showParticle, 
        pCoords = pData[[frame]]["Coords"];
        pField = Lookup[pData[[frame]]["Data"], renderTarget, ConstantArray[0., Length[pCoords]]];
        minV = Min[minV, Min[pField]]; maxV = Max[maxV, Max[pField]];
      ];
      If[hasGrid && showGrid,
        gCoords = gData[[frame]]["Coords"];
        gField = Lookup[gData[[frame]]["Data"], renderTarget, ConstantArray[0., Length[gCoords]]];
        minV = Min[minV, Min[gField]]; maxV = Max[maxV, Max[gField]];
      ];
      
      If[minV == Infinity, minV = 0.; maxV = 1.*^-6];
      If[minV == maxV, maxV = minV + 1.*^-6];
      
      If[hasGrid && showGrid,
        AppendTo[plotElements, {PointSize[gPtSize], Point[gCoords, VertexColors -> mpmPlot`Private`getFastColors[gField, minV, maxV]]}]
      ];
      If[showParticle,
        AppendTo[plotElements, {PointSize[pPtSize], Point[pCoords, VertexColors -> mpmPlot`Private`getFastColors[pField, minV, maxV]]}]
      ];
      
      Legended[
        Graphics[plotElements, PlotRange -> All, Frame -> True, FrameLabel -> {"X (m)", "Y (m)"}, 
                 PlotLabel -> StringForm["2D MPM | Frame ``/`` | Var: ``", frame, numFrames, renderTarget], ImageSize -> imgSz],
        BarLegend[{"Rainbow", {minV, maxV}}, LegendLabel -> renderTarget, Ticks -> mpmPlot`Private`sciTicks[minV, maxV]]
      ]
    ],
    
    (* ==================== Line 1 ==================== *)
    Row[{
      Control[{{renderTarget, defVar, "Output Var"}, vars}], 
      Spacer[30],
      Control[{{frame, 1, "Frame"}, 1, numFrames, 1, Appearance -> "Labeled"}]
    }],
    
    (* ==================== Line 2 ==================== *)
    Row[{
      Button["Export GIF", 
        Module[{outPath, frames},
          outPath = If[exportDir == "", FileNameJoin[{Directory[], ToString[exportName] <> ".gif"}], FileNameJoin[{exportDir, ToString[exportName] <> ".gif"}]];
          PrintTemporary["\|01f3ac Exporting dynamic-scaled 2D GIF..."];
          frames = Table[
            Module[{pC, pF, gC, gF, mn=Infinity, mx=-Infinity, elems={}},
              If[showParticle, 
                pC = pData[[i]]["Coords"]; pF = Lookup[pData[[i]]["Data"], renderTarget, ConstantArray[0., Length[pC]]];
                mn = Min[mn, Min[pF]]; mx = Max[mx, Max[pF]]];
              If[hasGrid && showGrid,
                gC = gData[[i]]["Coords"]; gF = Lookup[gData[[i]]["Data"], renderTarget, ConstantArray[0., Length[gC]]];
                mn = Min[mn, Min[gF]]; mx = Max[mx, Max[gF]]];
              If[mn == Infinity, mn = 0.; mx = 1.*^-6]; If[mn == mx, mx = mn + 1.*^-6];
              
              If[hasGrid && showGrid, AppendTo[elems, {PointSize[gPtSize], Point[gC, VertexColors -> mpmPlot`Private`getFastColors[gF, mn, mx]]}]];
              If[showParticle, AppendTo[elems, {PointSize[pPtSize], Point[pC, VertexColors -> mpmPlot`Private`getFastColors[pF, mn, mx]]}]];
              
              (* \:2705 \:6bcf\:4e00\:5e27\:4f7f\:7528\:5c40\:90e8\:6781\:503c (mn, mx)\:ff0c\:4f46 Ticks \:4f7f\:7528\:56fa\:5b9a\:5bbd\:5ea6\:7684\:79d1\:5b66\:8ba1\:6570\:6cd5 *)
              Legended[Graphics[elems, PlotRange -> All, Frame -> True, ImageSize -> 600], 
                BarLegend[{"Rainbow", {mn, mx}}, LegendLabel -> renderTarget, Ticks -> mpmPlot`Private`sciTicks[mn, mx]]]
            ], {i, 1, numFrames}];
          Export[outPath, frames, "GIF", "DisplayDurations" -> 0.08, "AnimationRepetitions" -> Infinity]; Print["\:2705 GIF Saved: ", outPath];
        ], Method -> "Queued", ImageSize -> {90, 25}],
      Spacer[10],
      InputField[Dynamic[exportName], String, FieldSize -> 15],
      ".gif"
    }],
    
    (* ==================== Line 3 ==================== *)
    Row[{
      Control[{{showParticle, True, "Particle"}, {True, False}}], Spacer[5],
      Control[{{pPtSize, OptionValue["pointSize"], "Size"}, 0.001, 0.02, 0.001, ImageSize -> Small}], 
      Spacer[30],
      Control[{{showGrid, True, "Grid"}, {True, False}, ControlType -> If[hasGrid, Checkbox, None]}], Spacer[If[hasGrid, 5, 0]],
      Control[{{gPtSize, OptionValue["pointSize"] * 0.8, "Size"}, 0.001, 0.02, 0.001, ImageSize -> Small, ControlType -> If[hasGrid, Manipulator, None]}]
    }],
    
    {{exportName, "2D_Animation"}, ControlType -> None},
    ControlPlacement -> Top,
    TrackedSymbols :> {frame, renderTarget, pPtSize, gPtSize, showParticle, showGrid}
  ]
];

mPlotData[pData_, opts : OptionsPattern[]] := mPlotData[pData, None, opts];

(* ========================================================================= *)
(* 5. 3D \:7ed8\:56fe\:6a21\:5757 (\:4e25\:683c\:4e09\:884c\:5168\:82f1\:6587\:7d27\:51d1\:5e03\:5c40 + \:52a8\:6001\:523b\:5ea6 + \:5b9a\:5bbd\:6392\:7248)               *)
(* ========================================================================= *)
Options[mPlotData3D] = {"defaultVar" -> "seq", "defaultView" -> "3D", "pointSize" -> 0.005, "exportPath" -> "", ColorFunction -> "Rainbow", ImageSize -> Large};

mPlotData3D[pData_, gData_, opts : OptionsPattern[]] /; !MatchQ[Unevaluated[gData], _Rule | _RuleDelayed] := Module[
  {hasGrid, numFrames, allHeaders, vars, defVar, defView, exportDir, imgSz},
  
  hasGrid = (Unevaluated[gData] =!= None);
  If[!MatchQ[Unevaluated[pData], _Symbol], MessageDialog["Please assign particle data to a variable first!"]; Return[$Failed]];
  If[hasGrid && !MatchQ[Unevaluated[gData], _Symbol], MessageDialog["Please assign grid data to a variable first!"]; Return[$Failed]];

  numFrames = If[hasGrid, Min[Length[pData], Length[gData]], Length[pData]];
  If[numFrames == 0, Return["Data is empty"]];
  
  allHeaders = DeleteDuplicates[Flatten[#["Headers"] & /@ pData]];
  If[hasGrid, allHeaders = DeleteDuplicates[Join[allHeaders, Flatten[#["Headers"] & /@ gData]]]];
  vars = Select[allHeaders, ! MemberQ[{"id", "type", "tag", "x", "y", "z"}, #] &];
  
  defVar = OptionValue["defaultVar"]; If[! MemberQ[vars, defVar], defVar = First[vars]];
  defView = OptionValue["defaultView"]; exportDir = OptionValue["exportPath"]; imgSz = OptionValue[ImageSize];

  Manipulate[
    Module[{pCoords, pField, gCoords, gField, minV, maxV, plotElements = {}, vParams},
      minV = Infinity; maxV = -Infinity;
      
      If[showParticle, 
        pCoords = pData[[frame]]["Coords"]; pField = Lookup[pData[[frame]]["Data"], renderTarget, ConstantArray[0., Length[pCoords]]];
        minV = Min[minV, Min[pField]]; maxV = Max[maxV, Max[pField]];
      ];
      If[hasGrid && showGrid,
        gCoords = gData[[frame]]["Coords"]; gField = Lookup[gData[[frame]]["Data"], renderTarget, ConstantArray[0., Length[gCoords]]];
        minV = Min[minV, Min[gField]]; maxV = Max[maxV, Max[gField]];
      ];
      
      If[minV == Infinity, minV = 0.; maxV = 1.*^-6]; If[minV == maxV, maxV = minV + 1.*^-6];
      
      If[hasGrid && showGrid,
        AppendTo[plotElements, {PointSize[gPtSize], Point[gCoords, VertexColors -> mpmPlot`Private`getFastColors[gField, minV, maxV]]}]
      ];
      If[showParticle,
        AppendTo[plotElements, {PointSize[pPtSize], Point[pCoords, VertexColors -> mpmPlot`Private`getFastColors[pField, minV, maxV]]}]
      ];
      
      vParams = Switch[viewOption,
        "3D", {{1.3, -2.4, 1.5}, {0, 0, 1}}, "XY", {{0, 0, Infinity}, {0, 1, 0}}, 
        "XZ", {{0, -Infinity, 0}, {0, 0, 1}}, "YZ", {{Infinity, 0, 0}, {0, 0, 1}}, 
        _, {{1.3, -2.4, 1.5}, {0, 0, 1}}
      ];

      Legended[
        Graphics3D[plotElements, PlotRange -> All, Axes -> True, AxesLabel -> {"X", "Y", "Z"}, Boxed -> True, 
                   ViewPoint -> vParams[[1]], ViewVertical -> vParams[[2]], 
                   PlotLabel -> StringForm["3D MPM | Frame ``/`` | Var: ``", frame, numFrames, renderTarget], ImageSize -> imgSz],
        BarLegend[{"Rainbow", {minV, maxV}}, LegendLabel -> renderTarget, Ticks -> mpmPlot`Private`sciTicks[minV, maxV]]
      ]
    ],
    
    (* ==================== Line 1 ==================== *)
    Row[{
      Control[{{renderTarget, defVar, "Output Var"}, vars}], 
      Spacer[30],
      Control[{{frame, 1, "Frame"}, 1, numFrames, 1, Appearance -> "Labeled"}]
    }],
    
    (* ==================== Line 2 ==================== *)
    Row[{
      Control[{{viewOption, defView, "View"}, {"3D", "XY", "XZ", "YZ"}}], 
      Spacer[30],
      Button["Export GIF", 
        Module[{outPath, frames, vP},
          outPath = If[exportDir == "", FileNameJoin[{Directory[], ToString[exportName] <> ".gif"}], FileNameJoin[{exportDir, ToString[exportName] <> ".gif"}]];
          PrintTemporary["\|01f3ac Exporting dynamic-scaled 3D GIF..."];
          vP = Switch[viewOption, "3D", {{1.3, -2.4, 1.5}, {0, 0, 1}}, "XY", {{0, 0, Infinity}, {0, 1, 0}}, "XZ", {{0, -Infinity, 0}, {0, 0, 1}}, "YZ", {{Infinity, 0, 0}, {0, 0, 1}}, _, {{1.3, -2.4, 1.5}, {0, 0, 1}}];
          
          (* \:2705 \:6bcf\:4e00\:5e27\:4f7f\:7528\:5c40\:90e8\:6781\:503c (mn, mx) \:91cd\:65b0\:8ba1\:7b97\:7f29\:653e *)
          frames = Table[
            Module[{pC, pF, gC, gF, mn=Infinity, mx=-Infinity, elems={}},
              If[showParticle, 
                pC = pData[[i]]["Coords"]; pF = Lookup[pData[[i]]["Data"], renderTarget, ConstantArray[0., Length[pC]]];
                mn = Min[mn, Min[pF]]; mx = Max[mx, Max[pF]]];
              If[hasGrid && showGrid,
                gC = gData[[i]]["Coords"]; gF = Lookup[gData[[i]]["Data"], renderTarget, ConstantArray[0., Length[gC]]];
                mn = Min[mn, Min[gF]]; mx = Max[mx, Max[gF]]];
              If[mn == Infinity, mn = 0.; mx = 1.*^-6]; If[mn == mx, mx = mn + 1.*^-6];
              
              If[hasGrid && showGrid, AppendTo[elems, {PointSize[gPtSize], Point[gC, VertexColors -> mpmPlot`Private`getFastColors[gF, mn, mx]]}]];
              If[showParticle, AppendTo[elems, {PointSize[pPtSize], Point[pC, VertexColors -> mpmPlot`Private`getFastColors[pF, mn, mx]]}]];
              
              Legended[Graphics3D[elems, PlotRange -> All, Axes -> True, Boxed -> True, ViewPoint -> vP[[1]], ViewVertical -> vP[[2]], ImageSize -> 600], 
                BarLegend[{"Rainbow", {mn, mx}}, LegendLabel -> renderTarget, Ticks -> mpmPlot`Private`sciTicks[mn, mx]]]
            ], {i, 1, numFrames}];
          Export[outPath, frames, "GIF", "DisplayDurations" -> 0.08, "AnimationRepetitions" -> Infinity]; Print["\:2705 GIF Saved: ", outPath];
        ], Method -> "Queued", ImageSize -> {90, 25}],
      Spacer[10],
      InputField[Dynamic[exportName], String, FieldSize -> 15],
      ".gif"
    }],
    
    (* ==================== Line 3 ==================== *)
    Row[{
      Control[{{showParticle, True, "Particle"}, {True, False}}], Spacer[5],
      Control[{{pPtSize, OptionValue["pointSize"], "Size"}, 0.001, 0.02, 0.001, ImageSize -> Small}], 
      Spacer[30],
      Control[{{showGrid, True, "Grid"}, {True, False}, ControlType -> If[hasGrid, Checkbox, None]}], Spacer[If[hasGrid, 5, 0]],
      Control[{{gPtSize, OptionValue["pointSize"] * 0.8, "Size"}, 0.001, 0.02, 0.001, ImageSize -> Small, ControlType -> If[hasGrid, Manipulator, None]}]
    }],
    
    {{exportName, "3D_Animation"}, ControlType -> None},
    ControlPlacement -> Top,
    TrackedSymbols :> {frame, renderTarget, viewOption, pPtSize, gPtSize, showParticle, showGrid}
  ]
];

mPlotData3D[pData_, opts : OptionsPattern[]] := mPlotData3D[pData, None, opts];

End[];
EndPackage[];

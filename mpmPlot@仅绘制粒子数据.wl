(* ::Package:: *)

(* ::Package:: *)
(**)


BeginPackage["mpmPlot`"];

(* ==================== \:63a5\:53e3\:4e0e Usage \:8bf4\:660e ==================== *)
loadMPMData::usage = 
  "loadMPMData[dirPath_String]\n" <>
  "\:3010\:8f93\:5165\:3011: dirPath (2D dump_p \:6587\:4ef6\:5939\:8def\:5f84)\:3002\n" <>
  "\:3010\:8f93\:51fa\:3011: 2D \:6570\:636e\:7ed3\:6784 Association \:5217\:8868\:3002";

mPlotData::usage = 
  "mPlotData[dataVar, opts___]\n" <>
  "\:3010\:8f93\:5165\:3011: dataVar (\:5fc5\:987b\:662f\:4fdd\:5b58\:4e86 2D \:6570\:636e\:7684\:53d8\:91cf\:540d\:ff0c\:5982 myData)\:3002\n" <>
  "\:3010\:9009\:9879\:3011: \"defaultVar\", \"pointSize\", ColorFunction, ImageSize, \"exportPath\"\:3002";

loadMPMData3D::usage = 
  "loadMPMData3D[dirPath_String]\n" <>
  "\:3010\:8f93\:5165\:3011: dirPath (3D dump_p \:6587\:4ef6\:5939\:8def\:5f84)\:3002\n" <>
  "\:3010\:8f93\:51fa\:3011: 3D \:6570\:636e\:7ed3\:6784 Association \:5217\:8868\:3002";

mPlotData3D::usage = 
  "mPlotData3D[dataVar, opts___]\n" <>
  "\:3010\:8f93\:5165\:3011: dataVar (\:5fc5\:987b\:662f\:4fdd\:5b58\:4e86 3D \:6570\:636e\:7684\:53d8\:91cf\:540d\:ff0c\:5982 myData)\:3002\n" <>
  "\:3010\:9009\:9879\:3011: \"defaultVar\", \"defaultView\", \"pointSize\", ColorFunction, ImageSize, \"exportPath\"\:3002";

Begin["`Private`"];

(* ========================================================================= *)
(* \:6838\:5fc3\:4fee\:590d\:ff1a\:5f3a\:5236\:4f20\:5740 (\:6309\:5f15\:7528\:4f20\:9012)                                             *)
(* \:963b\:6b62 Manipulate \:5c06\:6570\:767e\:5146\:6570\:636e\:5e8f\:5217\:5316\:8fdb Notebook UI \:4e2d\:ff0c\:6839\:9664\:5361\:987f\:5047\:6b7b\:73b0\:8c61           *)
(* ========================================================================= *)
SetAttributes[{mPlotData, mPlotData3D}, HoldFirst];

If[Kernels[] == {}, LaunchKernels[]];

(* ========================================================================= *)
(* 1. \:9ad8\:6027\:80fd\:989c\:8272\:5411\:91cf\:5316\:6620\:5c04 (\:767e\:4e07\:7ea7\:8d28\:70b9\:989c\:8272\:6781\:901f\:8ba1\:7b97)                              *)
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
(* 3. \:6570\:636e\:52a0\:8f7d\:51fd\:6570                                                           *)
(* ========================================================================= *)
loadMPMData[dataDir_String] := Module[{files, res},
  files = SortBy[FileNames["dump_p.*.LAMMPS", dataDir], ToExpression[StringExtract[FileNameTake[#], "." -> -2]] &];
  If[Length[files] == 0, MessageDialog["2D \:76ee\:5f55\:4e2d\:65e0 dump_p \:6587\:4ef6\:ff01"]; Return[$Failed]];
  PrintTemporary["\:26a1 \:6b63\:5728\:5e76\:884c\:52a0\:8f7d 2D \:6570\:636e..."];
  res = ParallelMap[fastParseRobust[#, False] &, files, DistributedContexts -> Automatic];
  res = Select[res, # =!= <||> &];
  Print["\:2705 2D \:6570\:636e\:5bfc\:5165\:6210\:529f\:ff0c\:5171 ", Length[res], " \:5e27\:ff01"];
  res
];

loadMPMData3D[dataDir_String] := Module[{files, res},
  files = SortBy[FileNames["dump_p.*.LAMMPS", dataDir], ToExpression[StringExtract[FileNameTake[#], "." -> -2]] &];
  If[Length[files] == 0, MessageDialog["3D \:76ee\:5f55\:4e2d\:65e0 dump_p \:6587\:4ef6\:ff01"]; Return[$Failed]];
  PrintTemporary["\:26a1 \:6b63\:5728\:5e76\:884c\:52a0\:8f7d 3D \:6570\:636e..."];
  res = ParallelMap[fastParseRobust[#, True] &, files, DistributedContexts -> Automatic];
  res = Select[res, # =!= <||> &];
  Print["\:2705 3D \:6570\:636e\:5bfc\:5165\:6210\:529f\:ff0c\:5171 ", Length[res], " \:5e27\:ff01"];
  res
];

(* ========================================================================= *)
(* 4. 2D \:7ed8\:56fe\:6a21\:5757 (\:4f20\:5740\:9632\:5361\:6b7b\:7248)                                               *)
(* ========================================================================= *)
Options[mPlotData] = {"defaultVar" -> "seq", "pointSize" -> 0.005, "exportPath" -> "", ColorFunction -> "Rainbow", ImageSize -> Large};

mPlotData[data_, opts : OptionsPattern[]] := Module[
  {vars, allHeaders, numFrames, defVar, exportDir, imgSz},
  
  (* \:8bed\:6cd5\:68c0\:67e5\:ff1a\:5f3a\:5236\:8981\:6c42\:4f20\:5165\:53d8\:91cf\:540d *)
  If[!MatchQ[Unevaluated[data], _Symbol],
    MessageDialog["\:274c \:6027\:80fd\:8b66\:544a\:ff1a\:4e3a\:4e86\:9632\:6b62\:524d\:7aef\:5361\:6b7b\:ff0c\:8bf7\:5148\:5c06\:6570\:636e\:8d4b\:503c\:7ed9\:53d8\:91cf\:ff01\n\n\:6b63\:786e\:7528\:6cd5\:ff1a\nmyData = loadMPMData[\"\:8def\:5f84\"];\nmPlotData[myData]"];
    Return[$Failed];
  ];
  
  numFrames = Length[data];
  If[numFrames == 0, Return["\:6570\:636e\:4e3a\:7a7a"]];
  
  allHeaders = DeleteDuplicates[Flatten[#["Headers"] & /@ data]];
  vars = Select[allHeaders, ! MemberQ[{"id", "type", "tag", "x", "y"}, #] &];
  defVar = OptionValue["defaultVar"]; If[! MemberQ[vars, defVar], defVar = First[vars]];
  exportDir = OptionValue["exportPath"]; imgSz = OptionValue[ImageSize];

  Manipulate[
    Module[{fObj, coords, fieldVal, minV, maxV, colors},
      fObj = data[[frame]]; (* \:56e0\:4e3a\:6709 HoldFirst\:ff0c\:8fd9\:91cc\:53ea\:4f1a\:8c03\:7528\:4f60\:7684\:5168\:5c40\:53d8\:91cf\:540d\:ff0c\:4e0d\:4f1a\:786c\:62f7\:8d1d *)
      coords = fObj["Coords"];
      fieldVal = Lookup[fObj["Data"], renderTarget, ConstantArray[0., Length[coords]]];
      minV = Min[fieldVal]; maxV = Max[fieldVal]; 
      colors = mpmPlot`Private`getFastColors[fieldVal, minV, maxV];
      
      Legended[
        Graphics[{PointSize[ptSize], Point[coords, VertexColors -> colors]}, 
          PlotRange -> All, Frame -> True, FrameLabel -> {"X (m)", "Y (m)"}, 
          PlotLabel -> StringForm["2D MPM | Frame ``/`` | Var: ``", frame, numFrames, renderTarget], ImageSize -> imgSz],
        BarLegend[{"Rainbow", {minV, If[minV == maxV, minV + 1.*^-6, maxV]}}, LegendLabel -> renderTarget]
      ]
    ],
    {{frame, 1, "\:5e27 (Frame)"}, 1, numFrames, 1, Appearance -> "Open"},
    {{renderTarget, defVar, "\:6e32\:67d3\:53d8\:91cf"}, vars},
    {{ptSize, OptionValue["pointSize"], "\:8d28\:70b9\:5927\:5c0f"}, 0.001, 0.02, 0.001},
    Item[Button["\|01f3ac \:5bfc\:51fa 2D GIF \:52a8\:753b", 
      Module[{out = If[exportDir == "", FileNameJoin[{Directory[], "2D_Animation.gif"}], exportDir], frames},
        PrintTemporary["\:6e32\:67d3 2D GIF \:4e2d..."];
        frames = Table[
          Module[{fv = Lookup[data[[i]]["Data"], renderTarget, ConstantArray[0., Length[data[[i]]["Coords"]]]], mn, mx, clrs},
            mn = Min[fv]; mx = Max[fv];
            clrs = mpmPlot`Private`getFastColors[fv, mn, mx];
            Legended[Graphics[{PointSize[ptSize], Point[data[[i]]["Coords"], VertexColors -> clrs]}, PlotRange -> All, Frame -> True, ImageSize -> 600], 
              BarLegend[{"Rainbow", {mn, If[mn == mx, mn + 1.*^-6, mx]}}]]
          ], {i, 1, numFrames}];
        Export[out, frames, "GIF", "DisplayDurations" -> 0.08, "AnimationRepetitions" -> Infinity];
        Print["\:2705 2D GIF \:6210\:529f\:5bfc\:51fa\:81f3: ", out];
      ], Method -> "Queued", ImageSize -> {200, 30}], Alignment -> Center],
    TrackedSymbols :> {frame, renderTarget, ptSize}
  ]
];

(* ========================================================================= *)
(* 5. 3D \:7ed8\:56fe\:6a21\:5757 (\:4f20\:5740\:9632\:5361\:6b7b\:7248)                                               *)
(* ========================================================================= *)
Options[mPlotData3D] = {"defaultVar" -> "seq", "defaultView" -> "3D \:89c6\:89d2", "pointSize" -> 0.005, "exportPath" -> "", ColorFunction -> "Rainbow", ImageSize -> Large};

mPlotData3D[data_, opts : OptionsPattern[]] := Module[
  {vars, allHeaders, numFrames, defVar, defView, exportDir, imgSz},
  
  If[!MatchQ[Unevaluated[data], _Symbol],
    MessageDialog["\:274c \:6027\:80fd\:8b66\:544a\:ff1a\:4e3a\:4e86\:9632\:6b62\:524d\:7aef\:5361\:6b7b\:ff0c\:8bf7\:5148\:5c06\:6570\:636e\:8d4b\:503c\:7ed9\:53d8\:91cf\:ff01\n\n\:6b63\:786e\:7528\:6cd5\:ff1a\nmyData = loadMPMData3D[\"\:8def\:5f84\"];\nmPlotData3D[myData]"];
    Return[$Failed];
  ];

  numFrames = Length[data];
  If[numFrames == 0, Return["\:6570\:636e\:4e3a\:7a7a"]];
  
  allHeaders = DeleteDuplicates[Flatten[#["Headers"] & /@ data]];
  vars = Select[allHeaders, ! MemberQ[{"id", "type", "tag", "x", "y", "z"}, #] &];
  defVar = OptionValue["defaultVar"]; If[! MemberQ[vars, defVar], defVar = First[vars]];
  defView = OptionValue["defaultView"]; exportDir = OptionValue["exportPath"]; imgSz = OptionValue[ImageSize];

  Manipulate[
    Module[{fObj, coords, fieldVal, minV, maxV, colors, vParams},
      fObj = data[[frame]];
      coords = fObj["Coords"];
      fieldVal = Lookup[fObj["Data"], renderTarget, ConstantArray[0., Length[coords]]];
      minV = Min[fieldVal]; maxV = Max[fieldVal]; 
      colors = mpmPlot`Private`getFastColors[fieldVal, minV, maxV];
      
      vParams = Switch[viewOption,
        "3D \:89c6\:89d2", {{1.3, -2.4, 1.5}, {0, 0, 1}}, "XY \:5e73\:9762 (\:4fef\:89c6)", {{0, 0, Infinity}, {0, 1, 0}},
        "XZ \:5e73\:9762 (\:6b63\:89c6)", {{0, -Infinity, 0}, {0, 0, 1}}, "YZ \:5e73\:9762 (\:4fa7\:89c6)", {{Infinity, 0, 0}, {0, 0, 1}}, _, {{1.3, -2.4, 1.5}, {0, 0, 1}}
      ];

      Legended[
        Graphics3D[{PointSize[ptSize], Point[coords, VertexColors -> colors]}, 
          PlotRange -> All, Axes -> True, AxesLabel -> {"X", "Y", "Z"}, Boxed -> True, 
          ViewPoint -> vParams[[1]], ViewVertical -> vParams[[2]], 
          PlotLabel -> StringForm["3D MPM | Frame ``/`` | Var: ``", frame, numFrames, renderTarget], ImageSize -> imgSz],
        BarLegend[{"Rainbow", {minV, If[minV == maxV, minV + 1.*^-6, maxV]}}, LegendLabel -> renderTarget]
      ]
    ],
    {{frame, 1, "\:5e27 (Frame)"}, 1, numFrames, 1, Appearance -> "Open"},
    {{renderTarget, defVar, "\:6e32\:67d3\:53d8\:91cf"}, vars},
    {{viewOption, defView, "\:89c6\:89d2\:9009\:62e9"}, {"3D \:89c6\:89d2", "XY \:5e73\:9762 (\:4fef\:89c6)", "XZ \:5e73\:9762 (\:6b63\:89c6)", "YZ \:5e73\:9762 (\:4fa7\:89c6)"}},
    {{ptSize, OptionValue["pointSize"], "\:8d28\:70b9\:5927\:5c0f"}, 0.001, 0.02, 0.001},
    Item[Button["\|01f3ac \:5bfc\:51fa 3D GIF \:52a8\:753b", 
      Module[{out = If[exportDir == "", FileNameJoin[{Directory[], "3D_Animation.gif"}], exportDir], frames, vP},
        PrintTemporary["\:6e32\:67d3 3D GIF \:4e2d..."];
        vP = Switch[viewOption,
          "3D \:89c6\:89d2", {{1.3, -2.4, 1.5}, {0, 0, 1}}, "XY \:5e73\:9762 (\:4fef\:89c6)", {{0, 0, Infinity}, {0, 1, 0}},
          "XZ \:5e73\:9762 (\:6b63\:89c6)", {{0, -Infinity, 0}, {0, 0, 1}}, "YZ \:5e73\:9762 (\:4fa7\:89c6)", {{Infinity, 0, 0}, {0, 0, 1}}, _, {{1.3, -2.4, 1.5}, {0, 0, 1}}
        ];
        frames = Table[
          Module[{fv = Lookup[data[[i]]["Data"], renderTarget, ConstantArray[0., Length[data[[i]]["Coords"]]]], mn, mx, clrs},
            mn = Min[fv]; mx = Max[fv]; 
            clrs = mpmPlot`Private`getFastColors[fv, mn, mx];
            Legended[Graphics3D[{PointSize[ptSize], Point[data[[i]]["Coords"], VertexColors -> clrs]}, 
              PlotRange -> All, Axes -> True, Boxed -> True, ViewPoint -> vP[[1]], ViewVertical -> vP[[2]], ImageSize -> 600], 
              BarLegend[{"Rainbow", {mn, If[mn == mx, mn + 1.*^-6, mx]}}]]
          ], {i, 1, numFrames}];
        Export[out, frames, "GIF", "DisplayDurations" -> 0.08, "AnimationRepetitions" -> Infinity];
        Print["\:2705 3D GIF \:6210\:529f\:5bfc\:51fa\:81f3: ", out];
      ], Method -> "Queued", ImageSize -> {200, 30}], Alignment -> Center],
    TrackedSymbols :> {frame, renderTarget, viewOption, ptSize}
  ]
];

End[];
EndPackage[];

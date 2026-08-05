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
(* 1. \:9ad8\:6027\:80fd\:989c\:8272\:5411\:91cf\:5316\:6620\:5c04                                                     *)
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
(* 4. 2D \:7ed8\:56fe\:6a21\:5757 (\:652f\:6301\:5355/\:53cc\:53d8\:91cf\:91cd\:8f7d + \:72ec\:7acb\:5c3a\:5bf8 + \:590d\:9009\:6846\:663e\:9690)                  *)
(* ========================================================================= *)
Options[mPlotData] = {"defaultVar" -> "seq", "pointSize" -> 0.005, "exportPath" -> "", ColorFunction -> "Rainbow", ImageSize -> Large};

mPlotData[pData_, gData_, opts : OptionsPattern[]] /; !MatchQ[Unevaluated[gData], _Rule | _RuleDelayed] := Module[
  {hasGrid, numFrames, allHeaders, vars, defVar, exportDir, imgSz},
  
  hasGrid = (Unevaluated[gData] =!= None);
  If[!MatchQ[Unevaluated[pData], _Symbol], MessageDialog["\:274c \:8bf7\:5c06\:7c92\:5b50\:6570\:636e\:8d4b\:503c\:7ed9\:53d8\:91cf\:ff01"]; Return[$Failed]];
  If[hasGrid && !MatchQ[Unevaluated[gData], _Symbol], MessageDialog["\:274c \:8bf7\:5c06\:7f51\:683c\:6570\:636e\:8d4b\:503c\:7ed9\:53d8\:91cf\:ff01"]; Return[$Failed]];
  
  numFrames = If[hasGrid, Min[Length[pData], Length[gData]], Length[pData]];
  If[numFrames == 0, Return["\:6570\:636e\:4e3a\:7a7a"]];
  
  allHeaders = DeleteDuplicates[Flatten[#["Headers"] & /@ pData]];
  If[hasGrid, allHeaders = DeleteDuplicates[Join[allHeaders, Flatten[#["Headers"] & /@ gData]]]];
  vars = Select[allHeaders, ! MemberQ[{"id", "type", "tag", "x", "y"}, #] &];
  
  defVar = OptionValue["defaultVar"]; If[! MemberQ[vars, defVar], defVar = First[vars]];
  exportDir = OptionValue["exportPath"]; imgSz = OptionValue[ImageSize];

  Manipulate[
    Module[{pCoords, pField, gCoords, gField, minV, maxV, plotElements = {}},
      minV = Infinity; maxV = -Infinity;
      
      (* \:6839\:636e\:590d\:9009\:6846\:72b6\:6001\:6309\:9700\:63d0\:53d6\:5e76\:8ba1\:7b97\:6781\:503c *)
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
      
      (* \:6839\:636e\:590d\:9009\:6846\:72b6\:6001\:53e0\:52a0\:56fe\:5c42 *)
      If[hasGrid && showGrid,
        AppendTo[plotElements, {PointSize[gPtSize], Point[gCoords, VertexColors -> mpmPlot`Private`getFastColors[gField, minV, maxV]]}]
      ];
      If[showParticle,
        AppendTo[plotElements, {PointSize[pPtSize], Point[pCoords, VertexColors -> mpmPlot`Private`getFastColors[pField, minV, maxV]]}]
      ];
      
      Legended[
        Graphics[plotElements, PlotRange -> All, Frame -> True, FrameLabel -> {"X (m)", "Y (m)"}, 
                 PlotLabel -> StringForm["2D MPM | Frame ``/`` | Var: ``", frame, numFrames, renderTarget], ImageSize -> imgSz],
        BarLegend[{"Rainbow", {minV, maxV}}, LegendLabel -> renderTarget]
      ]
    ],
    {{frame, 1, "\:5e27 (Frame)"}, 1, numFrames, 1, Appearance -> "Open"},
    {{renderTarget, defVar, "\:6e32\:67d3\:53d8\:91cf"}, vars},
    
    (* \:2705 \:590d\:9009\:6846\:663e\:793a\:63a7\:5236 *)
    {{showParticle, True, "\:663e\:793a\:7c92\:5b50 (Particle)"}, {True, False}},
    {{showGrid, True, "\:663e\:793a\:7f51\:683c (Grid)"}, {True, False}, ControlType -> If[hasGrid, Checkbox, None]},
     
    {{pPtSize, OptionValue["pointSize"], "\:7c92\:5b50\:5927\:5c0f (Particle)"}, 0.001, 0.02, 0.001},
    {{gPtSize, OptionValue["pointSize"] * 0.8, "\:7f51\:683c\:5927\:5c0f (Grid)"}, 0.001, 0.02, 0.001, ControlType -> If[hasGrid, Automatic, None]},
    
    Item[Button["\|01f3ac \:5bfc\:51fa 2D GIF \:52a8\:753b", 
      Module[{out = If[exportDir == "", FileNameJoin[{Directory[], "2D_Animation.gif"}], exportDir], frames},
        PrintTemporary["\:6e32\:67d3 2D GIF \:4e2d..."];
        frames = Table[
          Module[{pC, pF, gC, gF, mn=Infinity, mx=-Infinity, elems={}},
            If[showParticle, 
              pC = pData[[i]]["Coords"]; pF = Lookup[pData[[i]]["Data"], renderTarget, ConstantArray[0., Length[pC]]];
              mn = Min[mn, Min[pF]]; mx = Max[mx, Max[pF]];
            ];
            If[hasGrid && showGrid,
              gC = gData[[i]]["Coords"]; gF = Lookup[gData[[i]]["Data"], renderTarget, ConstantArray[0., Length[gC]]];
              mn = Min[mn, Min[gF]]; mx = Max[mx, Max[gF]];
            ];
            If[mn == Infinity, mn = 0.; mx = 1.*^-6]; If[mn == mx, mx = mn + 1.*^-6];
            If[hasGrid && showGrid, AppendTo[elems, {PointSize[gPtSize], Point[gC, VertexColors -> mpmPlot`Private`getFastColors[gF, mn, mx]]}]];
            If[showParticle, AppendTo[elems, {PointSize[pPtSize], Point[pC, VertexColors -> mpmPlot`Private`getFastColors[pF, mn, mx]]}]];
            Legended[Graphics[elems, PlotRange -> All, Frame -> True, ImageSize -> 600], BarLegend[{"Rainbow", {mn, mx}}]]
          ], {i, 1, numFrames}];
        Export[out, frames, "GIF", "DisplayDurations" -> 0.08, "AnimationRepetitions" -> Infinity]; Print["\:2705 2D GIF \:5bfc\:51fa\:81f3: ", out];
      ], Method -> "Queued", ImageSize -> {200, 30}], Alignment -> Center],
    TrackedSymbols :> {frame, renderTarget, pPtSize, gPtSize, showParticle, showGrid}
  ]
];

mPlotData[pData_, opts : OptionsPattern[]] := mPlotData[pData, None, opts];

(* ========================================================================= *)
(* 5. 3D \:7ed8\:56fe\:6a21\:5757 (\:652f\:6301\:5355/\:53cc\:53d8\:91cf\:91cd\:8f7d + \:72ec\:7acb\:5c3a\:5bf8 + \:590d\:9009\:6846\:663e\:9690)                  *)
(* ========================================================================= *)
Options[mPlotData3D] = {"defaultVar" -> "seq", "defaultView" -> "3D \:89c6\:89d2", "pointSize" -> 0.005, "exportPath" -> "", ColorFunction -> "Rainbow", ImageSize -> Large};

mPlotData3D[pData_, gData_, opts : OptionsPattern[]] /; !MatchQ[Unevaluated[gData], _Rule | _RuleDelayed] := Module[
  {hasGrid, numFrames, allHeaders, vars, defVar, defView, exportDir, imgSz},
  
  hasGrid = (Unevaluated[gData] =!= None);
  If[!MatchQ[Unevaluated[pData], _Symbol], MessageDialog["\:274c \:8bf7\:5c06\:7c92\:5b50\:6570\:636e\:8d4b\:503c\:7ed9\:53d8\:91cf\:ff01"]; Return[$Failed]];
  If[hasGrid && !MatchQ[Unevaluated[gData], _Symbol], MessageDialog["\:274c \:8bf7\:5c06\:7f51\:683c\:6570\:636e\:8d4b\:503c\:7ed9\:53d8\:91cf\:ff01"]; Return[$Failed]];

  numFrames = If[hasGrid, Min[Length[pData], Length[gData]], Length[pData]];
  If[numFrames == 0, Return["\:6570\:636e\:4e3a\:7a7a"]];
  
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
        "3D \:89c6\:89d2", {{1.3, -2.4, 1.5}, {0, 0, 1}}, "XY \:5e73\:9762 (\:4fef\:89c6)", {{0, 0, Infinity}, {0, 1, 0}},
        "XZ \:5e73\:9762 (\:6b63\:89c6)", {{0, -Infinity, 0}, {0, 0, 1}}, "YZ \:5e73\:9762 (\:4fa7\:89c6)", {{Infinity, 0, 0}, {0, 0, 1}}, _, {{1.3, -2.4, 1.5}, {0, 0, 1}}
      ];

      Legended[
        Graphics3D[plotElements, PlotRange -> All, Axes -> True, AxesLabel -> {"X", "Y", "Z"}, Boxed -> True, 
                   ViewPoint -> vParams[[1]], ViewVertical -> vParams[[2]], 
                   PlotLabel -> StringForm["3D MPM | Frame ``/`` | Var: ``", frame, numFrames, renderTarget], ImageSize -> imgSz],
        BarLegend[{"Rainbow", {minV, maxV}}, LegendLabel -> renderTarget]
      ]
    ],
    {{frame, 1, "\:5e27 (Frame)"}, 1, numFrames, 1, Appearance -> "Open"},
    {{renderTarget, defVar, "\:6e32\:67d3\:53d8\:91cf"}, vars},
    
    (* \:2705 \:590d\:9009\:6846\:663e\:793a\:63a7\:5236 *)
    {{showParticle, True, "\:663e\:793a\:7c92\:5b50 (Particle)"}, {True, False}},
    {{showGrid, True, "\:663e\:793a\:7f51\:683c (Grid)"}, {True, False}, ControlType -> If[hasGrid, Checkbox, None]},
    
    {{viewOption, defView, "\:89c6\:89d2\:9009\:62e9"}, {"3D \:89c6\:89d2", "XY \:5e73\:9762 (\:4fef\:89c6)", "XZ \:5e73\:9762 (\:6b63\:89c6)", "YZ \:5e73\:9762 (\:4fa7\:89c6)"}},
    
    {{pPtSize, OptionValue["pointSize"], "\:7c92\:5b50\:5927\:5c0f (Particle)"}, 0.001, 0.02, 0.001},
    {{gPtSize, OptionValue["pointSize"] * 0.8, "\:7f51\:683c\:5927\:5c0f (Grid)"}, 0.001, 0.02, 0.001, ControlType -> If[hasGrid, Automatic, None]},
    
    Item[Button["\|01f3ac \:5bfc\:51fa 3D GIF \:52a8\:753b", 
      Module[{out = If[exportDir == "", FileNameJoin[{Directory[], "3D_Animation.gif"}], exportDir], frames, vP},
        PrintTemporary["\:6e32\:67d3 3D GIF \:4e2d..."];
        vP = Switch[viewOption, "3D \:89c6\:89d2", {{1.3, -2.4, 1.5}, {0, 0, 1}}, "XY \:5e73\:9762 (\:4fef\:89c6)", {{0, 0, Infinity}, {0, 1, 0}}, "XZ \:5e73\:9762 (\:6b63\:89c6)", {{0, -Infinity, 0}, {0, 0, 1}}, "YZ \:5e73\:9762 (\:4fa7\:89c6)", {{Infinity, 0, 0}, {0, 0, 1}}, _, {{1.3, -2.4, 1.5}, {0, 0, 1}}];
        frames = Table[
          Module[{pC, pF, gC, gF, mn=Infinity, mx=-Infinity, elems={}},
            If[showParticle, 
              pC = pData[[i]]["Coords"]; pF = Lookup[pData[[i]]["Data"], renderTarget, ConstantArray[0., Length[pC]]];
              mn = Min[mn, Min[pF]]; mx = Max[mx, Max[pF]];
            ];
            If[hasGrid && showGrid,
              gC = gData[[i]]["Coords"]; gF = Lookup[gData[[i]]["Data"], renderTarget, ConstantArray[0., Length[gC]]];
              mn = Min[mn, Min[gF]]; mx = Max[mx, Max[gF]];
            ];
            If[mn == Infinity, mn = 0.; mx = 1.*^-6]; If[mn == mx, mx = mn + 1.*^-6];
            If[hasGrid && showGrid, AppendTo[elems, {PointSize[gPtSize], Point[gC, VertexColors -> mpmPlot`Private`getFastColors[gF, mn, mx]]}]];
            If[showParticle, AppendTo[elems, {PointSize[pPtSize], Point[pC, VertexColors -> mpmPlot`Private`getFastColors[pF, mn, mx]]}]];
            Legended[Graphics3D[elems, PlotRange -> All, Axes -> True, Boxed -> True, ViewPoint -> vP[[1]], ViewVertical -> vP[[2]], ImageSize -> 600], BarLegend[{"Rainbow", {mn, mx}}]]
          ], {i, 1, numFrames}];
        Export[out, frames, "GIF", "DisplayDurations" -> 0.08, "AnimationRepetitions" -> Infinity]; Print["\:2705 3D GIF \:5bfc\:51fa\:81f3: ", out];
      ], Method -> "Queued", ImageSize -> {200, 30}], Alignment -> Center],
    TrackedSymbols :> {frame, renderTarget, viewOption, pPtSize, gPtSize, showParticle, showGrid}
  ]
];

mPlotData3D[pData_, opts : OptionsPattern[]] := mPlotData3D[pData, None, opts];

End[];
EndPackage[];

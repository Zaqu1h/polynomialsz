(* Rodar/Run: wolframscript -file verif.wl *)

inputFile = "primesg5.txt";
outputFile = "resultado_" <> inputFile;

normalizeString[s_] := Module[{t = s},
  t = StringTrim[t];
  t = StringReplace[t, RegularExpression["\\s+"] -> ""];

  t = StringReplace[t, RegularExpression["iSqrt\\[([0-9]+)\\]"] -> "I*Sqrt[$1]"];
  t = StringReplace[t, RegularExpression["Sqrt\\[([0-9]+)\\]"] -> "Sqrt[$1]"];
  t = StringReplace[t, RegularExpression["Exp\\[\\-i\\*Pi/([0-9]+)\\]"] -> "Exp[-I*Pi/$1]"];
  t = StringReplace[t, RegularExpression["Exp\\[i\\*Pi/([0-9]+)\\]"] -> "Exp[I*Pi/$1]"];
  t = StringReplace[t, RegularExpression["Exp\\[(.*?)\\]"] -> "Exp[$1]"];
  t = StringReplace[t, "Pi" -> "\[Pi]"];

  t = StringReplace[t, RegularExpression["([0-9]+)x"] -> "$1*x"];
  t = StringReplace[t, RegularExpression["(^|[+\\-])x"] -> "$1 1*x"];
  t = StringReplace[t, RegularExpression["x(\\^[0-9]+)?\\("] -> "x$1*("];
  t = StringReplace[t, RegularExpression["([0-9]+)\\("] -> "$1*("];
  t = StringReplace[t, RegularExpression["\\)([0-9xI(])"] -> ")*$1"];
  t = StringReplace[t, ")(" -> ")*("];

  t = StringReplace[t, "x^1" -> "x"];
  t = StringReplace[t, "x^0" -> "1"];
  t = StringReplace[t, RegularExpression["([0-9.]+)I"] -> "$1*I"];
  t = StringReplace[t, "i" -> "I"];

  t
];

checkFactorization[lhsStr_, rhsStr_] := Catch[
  Module[{lhs, rhs, xs, tolerance = 10^-2, ratioRef = Indeterminate, 
          val1, val2, ratio},

    lhs = Quiet@Check[ToExpression[normalizeString[lhsStr]], $Failed];
    rhs = Quiet@Check[ToExpression[normalizeString[rhsStr]], $Failed];

    If[lhs === $Failed || rhs === $Failed, Throw[False]];

    xs = {-2, -1, 0, 1, 2, 7, 103, 0.5 + 0.5 I, -0.5 + 0.5 I};

    Do[
      val1 = N[lhs /. x -> xv, 30];
      val2 = N[rhs /. x -> xv, 30];

      If[Abs[val1] < tolerance && Abs[val2] < tolerance, Continue[]];
      If[(Abs[val1] < tolerance) != (Abs[val2] < tolerance), Throw[False]];

      ratio = val1/val2;

      If[ratioRef === Indeterminate,
        ratioRef = ratio,
        If[Abs[ratio - ratioRef] > 10^-3, Throw[False]]
      ];
    , {xv, xs}];

    True
  ]
];

processTXT[inputFile_, outputFile_] := Module[
  {lines, parts, poly, fact, status, correct = 0, incorrect = 0, stream},

  lines = ReadList[inputFile, String];
  stream = OpenWrite[outputFile];

  Do[
    If[StringTrim[line] == "", Continue[]];

    parts = StringSplit[line, "="];

    If[Length[parts] >= 2,
      poly = StringTrim[parts[[1]]];
      fact = StringTrim[parts[[-1]]];

        esc = FromCharacterCode[27];

        If[checkFactorization[poly, fact],
          status = "CORRETO";
          correct++;,
          status = "INCORRETO";
          incorrect++;
        ];

        WriteString[
          stream,
          poly, " = ", fact, " -> ", status, "\n"
        ];

        If[status === "INCORRETO",
          Print[
            poly, " = ", fact, " ->", esc, "[31m", status, esc, "[0m"
          ];
        ];

        (* Comentado para não mostrar os casos de acerto *)
        (*
        If[status === "CORRETO",
          Print[
            poly, " = ", fact, " -> ", esc, "[32m", status, esc, "[0m"
          ];
        ];
        *)
    ];

  , {line, lines}];

  WriteString[stream, "\n===== RESUMO =====\n"];
  WriteString[stream, "CORRETOS: ", correct, "\n"];
  WriteString[stream, "INCORRETOS: ", incorrect, "\n"];
  WriteString[stream, "TOTAL: ", correct + incorrect, "\n"];

  Close[stream];
  Print["Resultado salvo em: ", outputFile];
  Print["CORRETOS: ", correct];
  Print["INCORRETOS: ", incorrect];
  Print["TOTAL: ", correct + incorrect];
];

processTXT[inputFile, outputFile];

% Rodar/Run: octave verif.m

clear; clc;

inputFile = "primesg5.txt";
outputFile = ["resultado_" inputFile];

function t = normalizeString(s)

  t = strtrim(s);
  t = regexprep(t, '\s+', '');

  % Wolfram -> Octave
  t = regexprep(t, 'iSqrt\[([0-9]+)\]', '1i*sqrt($1)');
  t = regexprep(t, 'Sqrt\[([0-9]+)\]', 'sqrt($1)');
  t = regexprep(t, 'Exp\[\-i\*Pi/([0-9]+)\]', 'exp(-1i*pi/$1)');
  t = regexprep(t, 'Exp\[i\*Pi/([0-9]+)\]', 'exp(1i*pi/$1)');
  t = regexprep(t, 'Exp\[(.*?)\]', 'exp($1)');
  t = strrep(t, 'Pi', 'pi');

  % Multiplicacoes explicitas
  t = regexprep(t, '([0-9]+)x', '$1*x');
  t = regexprep(t, '(^|[+\-])x', '$1 1*x');
  t = regexprep(t, 'x(\^[0-9]+)?\(', 'x$1*(');
  t = regexprep(t, '([0-9]+)\(', '$1*(');
  t = regexprep(t, '\)([0-9x(])', ')*$1');
  t = strrep(t, ')(', ')*(');
  t = regexprep(t, '\)i', ')*i');

  % Simplificacoes
  t = strrep(t, 'x^1', 'x');
  t = strrep(t, 'x^0', '1');
  t = regexprep(t, '([0-9.]+)I', '$1i');

end

function ok = checkFactorization(lhs, rhs)

  tolerance = 1e-2;
  xs = [-2, -1, 0, 1, 2, 7, 103, 0.5+0.5i, -0.5+0.5i];

  ratio_ref = NaN;

  for k = 1:length(xs)
    try
      x = xs(k);
      val1 = eval(lhs);
      val2 = eval(rhs);

      % ignora ambos ~0
      if abs(val1) < tolerance && abs(val2) < tolerance
        continue;
      end

      % um zero e outro nao
      if (abs(val1) < tolerance) != (abs(val2) < tolerance)
        ok = false;
        return;
      end

      ratio = val1 / val2;

      if isnan(ratio_ref)
        ratio_ref = ratio;
      else
        if abs(ratio - ratio_ref) > 1e-3
          ok = false;
          return;
        end
      end

    catch
      ok = false;
      return;
    end
  end

  ok = true;
end

fid = fopen(inputFile, 'r');
if fid == -1
  error("Erro ao abrir arquivo.");
end

out = fopen(outputFile, 'w');

corretos = 0;
incorretos = 0;

fprintf("Processando arquivo...\n\n");

while ~feof(fid)
  line = strtrim(fgetl(fid));
  if isempty(line) || isempty(strfind(line, "="))
    continue;
  end

  parts = strsplit(line, '=');
  if length(parts) < 2
    continue;
  end

  lhsRaw = strtrim(parts{1});
  rhsRaw = strtrim(parts{end});

  lhs = normalizeString(lhsRaw);
  rhs = normalizeString(rhsRaw);

  if checkFactorization(lhs, rhs)
    status_txt = "CORRETO";
    status_term = "\033[32mCORRETO\033[0m";
    corretos++;

    fprintf(out, "%s = %s -> %s\n", lhsRaw, rhsRaw, status_txt);

    % Comentado para não mostrar os casos de acerto
    % fprintf("%s = %s -> %s\n", lhsRaw, rhsRaw, status_term);

  else
    status_txt = "INCORRETO";
    status_term = "\033[31mINCORRETO\033[0m";
    incorretos++;

    fprintf(out, "%s = %s -> %s\n", lhsRaw, rhsRaw, status_txt);

    fprintf("%s = %s -> %s\n", lhsRaw, rhsRaw, status_term);
  end
end

fclose(fid);

fprintf(out, "\n===== RESUMO =====\n");
fprintf(out, "CORRETOS: %d\n", corretos);
fprintf(out, "INCORRETOS: %d\n", incorretos);
fprintf(out, "TOTAL: %d\n", corretos + incorretos);

fclose(out);

fprintf("\n==== RESUMO ====\n");
fprintf("CORRETOS: %d\n", corretos);
fprintf("INCORRETOS: %d\n", incorretos);
fprintf("TOTAL: %d\n", corretos + incorretos);

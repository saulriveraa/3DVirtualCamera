function remart3D(h, s)

fn = h.FaceNormals;
to_del = sqrt(fn(:, :, 1).^2 + fn(:, :, 2).^2) > abs(fn(:, :, 3));
to_del = [to_del, false(s(1) - 1, 1); false(1, s(2))];
h.CData(to_del) = NaN;

end
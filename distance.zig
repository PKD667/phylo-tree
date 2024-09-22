// import stuuf
const std = @import("std");

fn compare(s: []const u8, b: []const u8) u32 {

    // m = len(s)
    const m = s.len;
    // n = len(b)
    const n = b.len;

    // v0 = [_]u32{0} ** m + 1
    var v0: [m + 1]u32 = undefined;

    // v1 = [_]u32{0} ** m + 1
    var v1: [m + 1]u32 = undefined;

    // for i in rnage 0..m
    for (s) |i| {

        // v1[0] = i+1
        v1[0] = i + 1;

        // for j in range 0..n
        for (b) |j| {

            // cost = if s[i] == b[j] 0 else 1
            const cost = if (s[i] == b[j]) 0 else 1;

            // v1[j+1] = min(v1[j] + 1, v0[j+1] + 1, v0[j] + cost)
            v1[j + 1] = std.math.min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost);
        }

        // v0 = v1.copy()
        v0 = v1;
    }

    return v1[n];
}

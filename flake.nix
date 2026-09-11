{
  description = "Kepler Formal equivalence checker";

  inputs = {
    self.submodules = true;
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
      kepler-formal = pkgs.callPackage ./nix/package.nix { src = self; };
    in
    {
      packages.${system} = {
        inherit kepler-formal;
        default = kepler-formal;
      };

      apps.${system}.default = {
        type = "app";
        program = "${kepler-formal}/bin/kepler-formal";
        meta.description = "Kepler Formal CLI";
      };

      checks.${system}.installed = pkgs.runCommand "kepler-formal-installed-check"
        { nativeBuildInputs = [ pkgs.bash pkgs.coreutils pkgs.gnugrep ]; }
        ''
          bash ${./nix/check.sh} ${kepler-formal}
          touch "$out"
        '';
    };
}

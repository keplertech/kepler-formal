{
  description = "Kepler Formal equivalence checker";

  inputs = {
    self.submodules = true;
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      forAllSystems = nixpkgs.lib.genAttrs [ "x86_64-linux" "aarch64-darwin" ];
    in
    {
      packages = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; };
          kepler-formal = pkgs.callPackage ./nix/package.nix { src = self; };
        in
        {
          inherit kepler-formal;
          default = kepler-formal;
        });

      apps = forAllSystems (system: {
        default = {
          type = "app";
          program = "${self.packages.${system}.kepler-formal}/bin/kepler-formal";
          meta.description = "Kepler Formal CLI";
        };
      });

      checks = forAllSystems (system:
        let pkgs = import nixpkgs { inherit system; };
        in {
          installed = pkgs.runCommand "kepler-formal-installed-check"
            { nativeBuildInputs = [ pkgs.bash pkgs.coreutils pkgs.gnugrep ]; }
            ''
              bash ${./nix/check.sh} ${self.packages.${system}.kepler-formal}
              touch "$out"
            '';
        });
    };
}

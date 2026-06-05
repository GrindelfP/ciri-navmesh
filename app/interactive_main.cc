/**
 * @file interactive_main.cc
 * @brief Interactive SDL2 viewer for navmesh triangulation and pathfinding.
 *
 * ## Controls
 *   - Type a filename in the text field and press Enter (or click Load)
 *   - Click DELA / GRDY / QGRD to triangulate with chosen algorithm
 *   - Left-click on the canvas: first click sets Start, second click sets Goal
 *     → path is computed and drawn automatically
 *   - Right-click: clear start/goal
 *   - Scroll wheel: zoom
 *   - Middle mouse drag: pan
 */

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>

#include "../lib/geometry/primitives.hpp"
#include "../lib/geometry/dcel.hh"
#include "../lib/triangulation/delaunay.hh"
#include "../lib/triangulation/greedy.hh"
#include "../lib/triangulation/quasi_greedy.hh"
#include "../lib/navmesh/navmesh.hh"
#include "../lib/navmesh/astar.hh"
#include "../lib/navmesh/funnel.hh"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
//  Constants / helpers
// ─────────────────────────────────────────────────────────────────────────────

static constexpr int WIN_W = 1100;
static constexpr int WIN_H = 800;
static constexpr int UI_H  = 60;   // top panel height
static constexpr int INFO_W = 220; // right info panel width

// Palette per algorithm (R,G,B)
static constexpr SDL_Color COL_DELA_EDGE  = {100, 140, 220, 255};
static constexpr SDL_Color COL_DELA_FACE  = { 30,  80, 160,  30};
static constexpr SDL_Color COL_GRDY_EDGE  = { 60, 160,  80, 255};
static constexpr SDL_Color COL_GRDY_FACE  = { 30, 140,  60,  30};
static constexpr SDL_Color COL_QGRD_EDGE  = {200, 120,  20, 255};
static constexpr SDL_Color COL_QGRD_FACE  = {180, 100,   0,  30};
static constexpr SDL_Color COL_PATH       = {220,  40,  40, 255};
static constexpr SDL_Color COL_CORRIDOR   = {255, 220,  40,  60};
static constexpr SDL_Color COL_START      = { 30, 180,  60, 255};
static constexpr SDL_Color COL_GOAL       = {200,  30,  30, 255};
static constexpr SDL_Color COL_POINT      = { 60,  60,  90, 255};
static constexpr SDL_Color COL_BG         = {245, 248, 252, 255};
static constexpr SDL_Color COL_UI_BG      = {220, 228, 238, 255};
static constexpr SDL_Color COL_BTN_NORM   = {180, 195, 215, 255};
static constexpr SDL_Color COL_BTN_HOV    = {140, 165, 200, 255};
static constexpr SDL_Color COL_BTN_ACT    = { 80, 130, 200, 255};
static constexpr SDL_Color COL_TEXT       = { 20,  20,  40, 255};
static constexpr SDL_Color COL_TEXT_LIGHT = {255, 255, 255, 255};
static constexpr SDL_Color COL_ERROR      = {180,  30,  30, 255};

// ─────────────────────────────────────────────────────────────────────────────
//  Tiny drawing helpers (no external deps)
// ─────────────────────────────────────────────────────────────────────────────

static void setColor(SDL_Renderer* r, SDL_Color c) {
    SDL_SetRenderDrawColor(r, c.r, c.g, c.b, c.a);
}

static void fillRect(SDL_Renderer* r, int x, int y, int w, int h, SDL_Color c) {
    setColor(r, c);
    SDL_Rect rect{x, y, w, h};
    SDL_RenderFillRect(r, &rect);
}

static void drawRect(SDL_Renderer* r, int x, int y, int w, int h, SDL_Color c) {
    setColor(r, c);
    SDL_Rect rect{x, y, w, h};
    SDL_RenderDrawRect(r, &rect);
}

// Filled circle via scanlines
static void fillCircle(SDL_Renderer* r, int cx, int cy, int rad, SDL_Color c) {
    setColor(r, c);
    for (int dy = -rad; dy <= rad; ++dy) {
        int dx = static_cast<int>(std::sqrt(rad * rad - dy * dy));
        SDL_RenderDrawLine(r, cx - dx, cy + dy, cx + dx, cy + dy);
    }
}

// Thick line (2px)
static void drawLine2(SDL_Renderer* r, int x0, int y0, int x1, int y1, SDL_Color c) {
    setColor(r, c);
    SDL_RenderDrawLine(r, x0,   y0,   x1,   y1);
    SDL_RenderDrawLine(r, x0+1, y0,   x1+1, y1);
    SDL_RenderDrawLine(r, x0,   y0+1, x1,   y1+1);
}

// Filled triangle (SDL has no native polygon fill — we use horizontal scan)
static void fillTriangle(SDL_Renderer* r,
                          int x0,int y0, int x1,int y1, int x2,int y2,
                          SDL_Color c)
{
    // Sort vertices by y
    if (y0 > y1) { std::swap(x0,x1); std::swap(y0,y1); }
    if (y0 > y2) { std::swap(x0,x2); std::swap(y0,y2); }
    if (y1 > y2) { std::swap(x1,x2); std::swap(y1,y2); }

    setColor(r, c);

    auto interpX = [](int xa, int ya, int xb, int yb, int y) -> int {
        if (ya == yb) return xa;
        return xa + (xb - xa) * (y - ya) / (yb - ya);
    };

    for (int y = y0; y <= y2; ++y) {
        int xA = interpX(x0,y0, x2,y2, y);
        int xB = (y < y1) ? interpX(x0,y0, x1,y1, y)
                           : interpX(x1,y1, x2,y2, y);
        if (xA > xB) std::swap(xA, xB);
        SDL_RenderDrawLine(r, xA, y, xB, y);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  TTF text helper (renders a string to a texture)
// ─────────────────────────────────────────────────────────────────────────────

struct TextCache {
    SDL_Renderer* ren{nullptr};
    TTF_Font*     font{nullptr};

    // Draw text at (x,y). Returns rendered width.
    int draw(const std::string& s, int x, int y, SDL_Color c) {
        if (!font) return 0;
        SDL_Surface* surf = TTF_RenderUTF8_Blended(font, s.c_str(), c);
        if (!surf) return 0;
        SDL_Texture* tex  = SDL_CreateTextureFromSurface(ren, surf);
        int w = surf->w, h = surf->h;
        SDL_FreeSurface(surf);
        if (!tex) return 0;
        SDL_Rect dst{x, y, w, h};
        SDL_RenderCopy(ren, tex, nullptr, &dst);
        SDL_DestroyTexture(tex);
        return w;
    }

    int textWidth(const std::string& s) {
        if (!font) return static_cast<int>(s.size()) * 8;
        int w = 0, h = 0;
        TTF_SizeUTF8(font, s.c_str(), &w, &h);
        return w;
    }
    int textHeight() {
        if (!font) return 14;
        return TTF_FontHeight(font);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  Viewport / coordinate transform
// ─────────────────────────────────────────────────────────────────────────────

struct Viewport {
    // Canvas area (pixels): [canvasX .. canvasX+canvasW] × [UI_H .. WIN_H]
    int canvasX{0}, canvasY{UI_H};
    int canvasW{WIN_W - INFO_W}, canvasH{WIN_H - UI_H};

    // World bounds
    double worldMinX{0}, worldMinY{0}, worldMaxX{1}, worldMaxY{1};

    // Pan / zoom (world coords of canvas centre)
    double cx{0.5}, cy{0.5};
    double scale{1.0}; // pixels per world unit

    void fitToWorld() {
        double wx = worldMaxX - worldMinX;
        double wy = worldMaxY - worldMinY;
        if (wx < 1e-9 || wy < 1e-9) return;
        double sx = canvasW / wx;
        double sy = canvasH / wy;
        scale = std::min(sx, sy) * 0.88;
        cx = (worldMinX + worldMaxX) * 0.5;
        cy = (worldMinY + worldMaxY) * 0.5;
    }

    // World → screen
    int toSX(double wx) const {
        return canvasX + canvasW/2 + static_cast<int>((wx - cx) * scale);
    }
    int toSY(double wy) const {
        // Y-flip: world Y up → screen Y down
        return canvasY + canvasH/2 - static_cast<int>((wy - cy) * scale);
    }

    // Screen → world
    double toWX(int sx) const {
        return cx + (sx - canvasX - canvasW/2.0) / scale;
    }
    double toWY(int sy) const {
        return cy - (sy - canvasY - canvasH/2.0) / scale;
    }

    bool inCanvas(int sx, int sy) const {
        return sx >= canvasX && sx < canvasX + canvasW
            && sy >= canvasY && sy < canvasY + canvasH;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  Button widget
// ─────────────────────────────────────────────────────────────────────────────

struct Button {
    SDL_Rect rect{};
    std::string label;
    bool active{false};   // currently selected / pressed

    bool contains(int x, int y) const {
        return x >= rect.x && x < rect.x+rect.w
            && y >= rect.y && y < rect.y+rect.h;
    }

    void draw(SDL_Renderer* r, TextCache& tc, bool hovered) const {
        SDL_Color bg = active   ? COL_BTN_ACT
                     : hovered  ? COL_BTN_HOV
                                : COL_BTN_NORM;
        fillRect(r, rect.x, rect.y, rect.w, rect.h, bg);
        drawRect(r, rect.x, rect.y, rect.w, rect.h, {100,120,150,255});

        SDL_Color tc_col = active ? COL_TEXT_LIGHT : COL_TEXT;
        int tw = tc.textWidth(label);
        int th = tc.textHeight();
        tc.draw(label,
                rect.x + (rect.w - tw)/2,
                rect.y + (rect.h - th)/2,
                tc_col);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  TextInput widget
// ─────────────────────────────────────────────────────────────────────────────

struct TextInput {
    SDL_Rect rect{};
    std::string text;
    bool focused{false};

    bool contains(int x, int y) const {
        return x >= rect.x && x < rect.x+rect.w
            && y >= rect.y && y < rect.y+rect.h;
    }

    void draw(SDL_Renderer* r, TextCache& tc) const {
        SDL_Color bg = focused ? SDL_Color{255,255,255,255}
                               : SDL_Color{235,240,248,255};
        fillRect(r, rect.x, rect.y, rect.w, rect.h, bg);
        drawRect(r, rect.x, rect.y, rect.w, rect.h,
                 focused ? SDL_Color{80,130,200,255}
                         : SDL_Color{160,170,185,255});

        std::string disp = text;
        if (focused) disp += "|";
        tc.draw(disp, rect.x+6, rect.y + (rect.h - tc.textHeight())/2, COL_TEXT);
    }

    void handleKey(SDL_Keycode key, const std::string& input_text) {
        if (!focused) return;
        if (key == SDLK_BACKSPACE && !text.empty())
            text.pop_back();
    }

    void appendText(const std::string& t) {
        if (focused) text += t;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  Application state
// ─────────────────────────────────────────────────────────────────────────────

enum class AlgoChoice { None, Delaunay, Greedy, Quasi };

struct AppState {
    // Data
    std::vector<geometry::Point2D> points;
    geometry::DCEL                 dcel;
    navmesh::NavMesh               nm;
    AlgoChoice                     algo{AlgoChoice::None};

    // Pathfinding
    std::optional<geometry::Point2D> startPt;
    std::optional<geometry::Point2D> goalPt;
    std::vector<navmesh::NodeIdx>    corridor;
    std::vector<geometry::Point2D>   waypoints;
    double pathLen{0.0};
    std::size_t turns{0};
    bool pathFound{false};

    // Stats
    std::size_t triCount{0};
    double      totalWeight{0.0};
    double      buildMs{0.0};
    std::size_t flipCount{0};

    // UI state
    std::string statusMsg;
    bool        statusIsError{false};

    void clearPath() {
        startPt.reset();
        goalPt.reset();
        corridor.clear();
        waypoints.clear();
        pathFound = false;
        pathLen = 0; turns = 0;
    }

    void clearTriangulation() {
        dcel.clear();
        nm.clear();
        triCount = 0; totalWeight = 0; buildMs = 0; flipCount = 0;
        clearPath();
        algo = AlgoChoice::None;
    }

    void runPath() {
        if (!startPt || !goalPt) return;
        corridor.clear();
        waypoints.clear();
        pathFound = false;

        navmesh::NodeIdx sn = nm.findNode(*startPt, dcel);
        navmesh::NodeIdx gn = nm.findNode(*goalPt,  dcel);

        if (sn == navmesh::kInvalidNode || gn == navmesh::kInvalidNode) {
            statusMsg = "Start or goal outside mesh!";
            statusIsError = true;
            return;
        }

        pathfinding::AStar astar;
        auto res = astar.search(nm, sn, gn);
        if (!res.found) {
            statusMsg = "No path found between selected points.";
            statusIsError = true;
            return;
        }

        corridor = res.nodes;
        pathfinding::FunnelAlgorithm funnel;
        auto smooth = funnel.smooth(*startPt, *goalPt, corridor, nm, dcel);
        waypoints = smooth.waypoints;
        pathLen   = smooth.totalLength;
        turns     = smooth.turnCount;
        pathFound = true;
        statusMsg = "Path found! Len=" + [&]{
            std::ostringstream ss;
            ss << std::fixed;
            ss.precision(2);
            ss << pathLen;
            return ss.str();
        }() + " Turns=" + std::to_string(turns);
        statusIsError = false;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
//  File reading (same logic as main.cc)
// ─────────────────────────────────────────────────────────────────────────────

static bool loadPoints(const std::string& path,
                        std::vector<geometry::Point2D>& out,
                        std::string& err)
{
    std::ifstream in(path);
    if (!in) { err = "Cannot open: " + path; return false; }

    out.clear();
    std::string line;
    int ln = 0;
    while (std::getline(in, line)) {
        ++ln;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream ss(line);
        double x{}, y{};
        if (!(ss >> x >> y)) {
            err = "Parse error line " + std::to_string(ln);
            return false;
        }
        out.push_back({x, y});
    }

    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());

    if (out.size() < 3) {
        err = "Need >= 3 distinct points, got " + std::to_string(out.size());
        return false;
    }
    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Triangulation runner
// ─────────────────────────────────────────────────────────────────────────────

static bool runTriangulation(AppState& state, AlgoChoice choice, std::string& err)
{
    state.clearTriangulation();

    try {
        triangulation::TriangulationResult res;

        if (choice == AlgoChoice::Delaunay) {
            triangulation::DelaunayTriangulator dt;
            res = dt.triangulate(state.points, state.dcel);
        } else if (choice == AlgoChoice::Greedy) {
            triangulation::GreedyTriangulator gt;
            res = gt.triangulate(state.points, state.dcel);
        } else {
            triangulation::QuasiGreedyTriangulator qg(/*doFlips=*/true);
            res = qg.triangulate(state.points, state.dcel);
            state.flipCount = res.flipCount;
        }

        state.nm.build(state.dcel);
        state.algo        = choice;
        state.triCount    = res.triangleCount;
        state.totalWeight = res.totalWeight;
        state.buildMs     = res.elapsed.count() * 1000.0;
        return true;
    } catch (const std::exception& e) {
        err = std::string("Triangulation failed: ") + e.what();
        return false;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Drawing: triangulation
// ─────────────────────────────────────────────────────────────────────────────

static void drawTriangulation(SDL_Renderer* r, const AppState& state,
                               const Viewport& vp)
{
    if (state.algo == AlgoChoice::None) return;

    SDL_Color faceCol, edgeCol;
    switch (state.algo) {
        case AlgoChoice::Delaunay:
            faceCol = COL_DELA_FACE; edgeCol = COL_DELA_EDGE; break;
        case AlgoChoice::Greedy:
            faceCol = COL_GRDY_FACE; edgeCol = COL_GRDY_EDGE; break;
        default:
            faceCol = COL_QGRD_FACE; edgeCol = COL_QGRD_EDGE; break;
    }

    const auto& dcel = state.dcel;

    // Enable blending for face fill
    SDL_SetRenderDrawBlendMode(r, SDL_BLENDMODE_BLEND);

    // Fill faces
    for (geometry::FaceIdx fi = 1; fi < dcel.faceCount(); ++fi) {
        if (dcel.face(fi).dead) continue;
        auto [va, vb, vc] = dcel.faceVertices(fi);
        int ax = vp.toSX(dcel.vertex(va).pos.x), ay = vp.toSY(dcel.vertex(va).pos.y);
        int bx = vp.toSX(dcel.vertex(vb).pos.x), by = vp.toSY(dcel.vertex(vb).pos.y);
        int cx = vp.toSX(dcel.vertex(vc).pos.x), cy = vp.toSY(dcel.vertex(vc).pos.y);
        fillTriangle(r, ax,ay, bx,by, cx,cy, faceCol);
    }

    SDL_SetRenderDrawBlendMode(r, SDL_BLENDMODE_NONE);

    // Draw edges
    for (geometry::HalfEdgeIdx hi = 0; hi < dcel.halfEdgeCount(); ++hi) {
        if (!dcel.isHalfEdgeLive(hi)) continue;
        const auto& he = dcel.halfEdge(hi);
        if (he.dead) continue;
        if (hi > he.twin) continue;
        if (he.face == geometry::kOuterFace &&
            dcel.halfEdge(he.twin).face == geometry::kOuterFace) continue;

        int x1 = vp.toSX(dcel.vertex(he.origin).pos.x);
        int y1 = vp.toSY(dcel.vertex(he.origin).pos.y);
        int x2 = vp.toSX(dcel.vertex(dcel.halfEdge(he.twin).origin).pos.x);
        int y2 = vp.toSY(dcel.vertex(dcel.halfEdge(he.twin).origin).pos.y);

        setColor(r, edgeCol);
        SDL_RenderDrawLine(r, x1, y1, x2, y2);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Drawing: corridor
// ─────────────────────────────────────────────────────────────────────────────

static void drawCorridor(SDL_Renderer* r, const AppState& state,
                          const Viewport& vp)
{
    if (state.corridor.empty()) return;
    SDL_SetRenderDrawBlendMode(r, SDL_BLENDMODE_BLEND);

    const auto& dcel = state.dcel;
    const auto& nm   = state.nm;

    for (navmesh::NodeIdx ni : state.corridor) {
        const auto& nd = nm.node(ni);
        auto [va, vb, vc] = dcel.faceVertices(nd.faceIdx);
        int ax = vp.toSX(dcel.vertex(va).pos.x), ay = vp.toSY(dcel.vertex(va).pos.y);
        int bx = vp.toSX(dcel.vertex(vb).pos.x), by = vp.toSY(dcel.vertex(vb).pos.y);
        int cx = vp.toSX(dcel.vertex(vc).pos.x), cy = vp.toSY(dcel.vertex(vc).pos.y);
        fillTriangle(r, ax,ay, bx,by, cx,cy, COL_CORRIDOR);
    }
    SDL_SetRenderDrawBlendMode(r, SDL_BLENDMODE_NONE);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Drawing: path
// ─────────────────────────────────────────────────────────────────────────────

static void drawPath(SDL_Renderer* r, const AppState& state, const Viewport& vp) {
    const auto& wp = state.waypoints;
    if (wp.size() < 2) return;

    for (std::size_t i = 0; i + 1 < wp.size(); ++i) {
        drawLine2(r,
                  vp.toSX(wp[i].x),   vp.toSY(wp[i].y),
                  vp.toSX(wp[i+1].x), vp.toSY(wp[i+1].y),
                  COL_PATH);
    }
    // Interior waypoint dots
    for (std::size_t i = 1; i + 1 < wp.size(); ++i) {
        fillCircle(r, vp.toSX(wp[i].x), vp.toSY(wp[i].y), 4, COL_PATH);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Drawing: points
// ─────────────────────────────────────────────────────────────────────────────

static void drawPoints(SDL_Renderer* r, const AppState& state, const Viewport& vp) {
    for (const auto& p : state.points) {
        fillCircle(r, vp.toSX(p.x), vp.toSY(p.y), 3, COL_POINT);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Drawing: start / goal markers
// ─────────────────────────────────────────────────────────────────────────────

static void drawMarker(SDL_Renderer* r, TextCache& tc,
                        int sx, int sy, SDL_Color col, const char* label)
{
    fillCircle(r, sx, sy, 10, col);
    // White ring
    setColor(r, {255,255,255,255});
    for (int a = 0; a < 360; a += 5) {
        double rad = a * M_PI / 180.0;
        int px = sx + static_cast<int>(10 * std::cos(rad));
        int py = sy + static_cast<int>(10 * std::sin(rad));
        SDL_RenderDrawPoint(r, px, py);
    }
    int tw = tc.textWidth(label);
    int th = tc.textHeight();
    tc.draw(label, sx - tw/2, sy - th/2, COL_TEXT_LIGHT);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Info panel (right side)
// ─────────────────────────────────────────────────────────────────────────────

static void drawInfoPanel(SDL_Renderer* r, TextCache& tc,
                           const AppState& state)
{
    int px = WIN_W - INFO_W + 8;
    int py = UI_H + 10;
    int lh = tc.textHeight() + 4;

    fillRect(r, WIN_W - INFO_W, UI_H, INFO_W, WIN_H - UI_H, {232,238,248,255});
    setColor(r, {190,200,215,255});
    SDL_RenderDrawLine(r, WIN_W - INFO_W, UI_H, WIN_W - INFO_W, WIN_H);

    auto line = [&](const std::string& s, SDL_Color c = COL_TEXT) {
        tc.draw(s, px, py, c);
        py += lh;
    };

    line("── Algorithm ──", {80,100,140,255});
    std::string algoStr = "None";
    if (state.algo == AlgoChoice::Delaunay) algoStr = "Delaunay (BW)";
    else if (state.algo == AlgoChoice::Greedy) algoStr = "Greedy MWT";
    else if (state.algo == AlgoChoice::Quasi) algoStr = "Quasi-Greedy";
    line(algoStr, COL_BTN_ACT);
    py += 4;

    if (state.algo != AlgoChoice::None) {
        line("Triangles: " + std::to_string(state.triCount));

        std::ostringstream wss;
        wss << std::fixed; wss.precision(2);
        wss << state.totalWeight;
        line("Weight:  " + wss.str());

        std::ostringstream bss;
        bss << std::fixed; bss.precision(1);
        bss << state.buildMs;
        line("Build:   " + bss.str() + " ms");

        if (state.algo == AlgoChoice::Quasi && state.flipCount > 0) {
            line("Flips:   " + std::to_string(state.flipCount));
        }
    }
    py += 8;
    line("── Path ──", {80,100,140,255});

    if (state.startPt) {
        std::ostringstream ss;
        ss << std::fixed; ss.precision(2);
        ss << "S(" << state.startPt->x << "," << state.startPt->y << ")";
        line(ss.str(), COL_START);
    } else {
        line("Start: (click)", {150,150,160,255});
    }

    if (state.goalPt) {
        std::ostringstream ss;
        ss << std::fixed; ss.precision(2);
        ss << "G(" << state.goalPt->x << "," << state.goalPt->y << ")";
        line(ss.str(), COL_GOAL);
    } else {
        line("Goal:  (click)", {150,150,160,255});
    }

    if (state.pathFound) {
        py += 4;
        std::ostringstream pss;
        pss << std::fixed; pss.precision(2);
        pss << state.pathLen;
        line("PathLen: " + pss.str());
        line("Turns:   " + std::to_string(state.turns));
    }

    py += 8;
    line("── Controls ──", {80,100,140,255});
    line("LClick: set S/G");
    line("RClick: clear");
    line("Scroll: zoom");
    line("MBtn+drag: pan");
}

// ─────────────────────────────────────────────────────────────────────────────
//  main
// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    (void)argc; (void)argv;

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cerr << "SDL_Init: " << SDL_GetError() << "\n";
        return 1;
    }
    if (TTF_Init() != 0) {
        std::cerr << "TTF_Init: " << TTF_GetError() << "\n";
        SDL_Quit();
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow(
        "NavMesh Interactive Viewer",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIN_W, WIN_H,
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    if (!window) {
        std::cerr << "CreateWindow: " << SDL_GetError() << "\n";
        TTF_Quit(); SDL_Quit();
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(
        window, -1,
        SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        std::cerr << "CreateRenderer: " << SDL_GetError() << "\n";
        SDL_DestroyWindow(window);
        TTF_Quit(); SDL_Quit();
        return 1;
    }
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    // ── Load font ─────────────────────────────────────────────────────────────
    // Try several common system font locations on macOS
    TTF_Font* font = nullptr;
    const std::vector<std::string> fontPaths = {
        "/System/Library/Fonts/Supplemental/Courier New.ttf",
        "/System/Library/Fonts/Monaco.ttf",
        "/System/Library/Fonts/Menlo.ttc",
        "/Library/Fonts/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/System/Library/Fonts/SFNSMono.ttf",
        "/System/Library/Fonts/SFNSDisplay.ttf",
    };
    for (const auto& fp : fontPaths) {
        font = TTF_OpenFont(fp.c_str(), 13);
        if (font) { std::cout << "Loaded font: " << fp << "\n"; break; }
    }
    if (!font) {
        std::cerr << "Warning: could not load any system font. "
                     "Text rendering disabled.\n";
    }

    TextCache tc{renderer, font};

    // ── UI layout ─────────────────────────────────────────────────────────────
    TextInput filenameInput;
    filenameInput.rect = {10, 10, 280, 36};
    filenameInput.text = "points.txt";
    filenameInput.focused = true;

    Button btnLoad  = {{300, 10, 70, 36}, "Load"};
    Button btnDela  = {{385, 10, 70, 36}, "DELA"};
    Button btnGrdy  = {{460, 10, 70, 36}, "GRDY"};
    Button btnQgrd  = {{535, 10, 70, 36}, "QGRD"};
    Button btnClear = {{620, 10, 80, 36}, "ClearPath"};

    AppState state;
    Viewport vp;
    vp.canvasW = WIN_W - INFO_W;

    std::string hoveredBtn = "";

    bool panActive = false;
    int panLastX = 0, panLastY = 0;

    bool running = true;
    SDL_Event ev;

    SDL_StartTextInput();

    while (running) {
        // ── Events ────────────────────────────────────────────────────────────
        while (SDL_PollEvent(&ev)) {
            switch (ev.type) {
            case SDL_QUIT:
                running = false;
                break;

            case SDL_TEXTINPUT:
                filenameInput.appendText(ev.text.text);
                break;

            case SDL_KEYDOWN: {
                SDL_Keycode key = ev.key.keysym.sym;
                if (key == SDLK_BACKSPACE && filenameInput.focused) {
                    if (!filenameInput.text.empty())
                        filenameInput.text.pop_back();
                } else if (key == SDLK_RETURN || key == SDLK_KP_ENTER) {
                    // Load file
                    std::string err;
                    if (loadPoints(filenameInput.text, state.points, err)) {
                        state.clearTriangulation();
                        // Set world bounds
                        double mnx = state.points[0].x, mxx = state.points[0].x;
                        double mny = state.points[0].y, mxy = state.points[0].y;
                        for (const auto& p : state.points) {
                            mnx = std::min(mnx, p.x); mxx = std::max(mxx, p.x);
                            mny = std::min(mny, p.y); mxy = std::max(mxy, p.y);
                        }
                        double pad = std::max(mxx-mnx, mxy-mny) * 0.1 + 1.0;
                        vp.worldMinX = mnx - pad; vp.worldMaxX = mxx + pad;
                        vp.worldMinY = mny - pad; vp.worldMaxY = mxy + pad;
                        vp.fitToWorld();
                        state.statusMsg = "Loaded " +
                            std::to_string(state.points.size()) + " points.";
                        state.statusIsError = false;
                    } else {
                        state.statusMsg = err;
                        state.statusIsError = true;
                    }
                } else if (key == SDLK_ESCAPE) {
                    filenameInput.focused = false;
                }
                break;
            }

            case SDL_MOUSEBUTTONDOWN: {
                int mx = ev.button.x, my = ev.button.y;

                // UI clicks
                if (ev.button.button == SDL_BUTTON_LEFT) {
                    filenameInput.focused = filenameInput.contains(mx, my);

                    if (btnLoad.contains(mx, my)) {
                        std::string err;
                        if (loadPoints(filenameInput.text, state.points, err)) {
                            state.clearTriangulation();
                            double mnx = state.points[0].x, mxx = state.points[0].x;
                            double mny = state.points[0].y, mxy = state.points[0].y;
                            for (const auto& p : state.points) {
                                mnx=std::min(mnx,p.x); mxx=std::max(mxx,p.x);
                                mny=std::min(mny,p.y); mxy=std::max(mxy,p.y);
                            }
                            double pad = std::max(mxx-mnx,mxy-mny)*0.1+1.0;
                            vp.worldMinX=mnx-pad; vp.worldMaxX=mxx+pad;
                            vp.worldMinY=mny-pad; vp.worldMaxY=mxy+pad;
                            vp.fitToWorld();
                            state.statusMsg = "Loaded " +
                                std::to_string(state.points.size()) + " points.";
                            state.statusIsError = false;
                        } else {
                            state.statusMsg = err;
                            state.statusIsError = true;
                        }
                    }
                    else if (btnDela.contains(mx, my)) {
                        if (!state.points.empty()) {
                            std::string err;
                            if (!runTriangulation(state, AlgoChoice::Delaunay, err)) {
                                state.statusMsg = err; state.statusIsError = true;
                            } else {
                                state.statusMsg = "Delaunay done.";
                                state.statusIsError = false;
                            }
                        }
                    }
                    else if (btnGrdy.contains(mx, my)) {
                        if (!state.points.empty()) {
                            std::string err;
                            if (!runTriangulation(state, AlgoChoice::Greedy, err)) {
                                state.statusMsg = err; state.statusIsError = true;
                            } else {
                                state.statusMsg = "Greedy MWT done.";
                                state.statusIsError = false;
                            }
                        }
                    }
                    else if (btnQgrd.contains(mx, my)) {
                        if (!state.points.empty()) {
                            std::string err;
                            if (!runTriangulation(state, AlgoChoice::Quasi, err)) {
                                state.statusMsg = err; state.statusIsError = true;
                            } else {
                                state.statusMsg = "Quasi-Greedy done.";
                                state.statusIsError = false;
                            }
                        }
                    }
                    else if (btnClear.contains(mx, my)) {
                        state.clearPath();
                        state.statusMsg = "Path cleared.";
                        state.statusIsError = false;
                    }
                    else if (vp.inCanvas(mx, my) &&
                             state.algo != AlgoChoice::None)
                    {
                        // Canvas click: set start or goal
                        geometry::Point2D wp{ vp.toWX(mx), vp.toWY(my) };
                        if (!state.startPt) {
                            state.startPt = wp;
                            state.statusMsg = "Start set. Now click goal.";
                            state.statusIsError = false;
                        } else if (!state.goalPt) {
                            state.goalPt = wp;
                            state.runPath();
                        } else {
                            // Both already set: restart from new start
                            state.clearPath();
                            state.startPt = wp;
                            state.statusMsg = "Start reset. Now click goal.";
                            state.statusIsError = false;
                        }
                    }
                }
                else if (ev.button.button == SDL_BUTTON_RIGHT) {
                    state.clearPath();
                    state.statusMsg = "Path cleared.";
                    state.statusIsError = false;
                }
                else if (ev.button.button == SDL_BUTTON_MIDDLE) {
                    panActive = true;
                    panLastX = mx; panLastY = my;
                }
                break;
            }

            case SDL_MOUSEBUTTONUP:
                if (ev.button.button == SDL_BUTTON_MIDDLE)
                    panActive = false;
                break;

            case SDL_MOUSEMOTION: {
                if (panActive && vp.inCanvas(ev.motion.x, ev.motion.y)) {
                    int dx = ev.motion.x - panLastX;
                    int dy = ev.motion.y - panLastY;
                    vp.cx -= dx / vp.scale;
                    vp.cy += dy / vp.scale;
                    panLastX = ev.motion.x;
                    panLastY = ev.motion.y;
                }
                break;
            }

            case SDL_MOUSEWHEEL: {
                int mx, my;
                SDL_GetMouseState(&mx, &my);
                if (vp.inCanvas(mx, my)) {
                    double factor = (ev.wheel.y > 0) ? 1.15 : (1.0/1.15);
                    // Zoom towards mouse position
                    double wx = vp.toWX(mx), wy = vp.toWY(my);
                    vp.scale *= factor;
                    // Adjust centre so the point under cursor stays fixed
                    vp.cx = wx - (mx - vp.canvasX - vp.canvasW/2.0) / vp.scale;
                    vp.cy = wy + (my - vp.canvasY - vp.canvasH/2.0) / vp.scale;
                }
                break;
            }

            case SDL_WINDOWEVENT:
                if (ev.window.event == SDL_WINDOWEVENT_RESIZED) {
                    // Re-query actual window size
                    SDL_GetWindowSize(window,
                        &vp.canvasW, &vp.canvasH);
                    // subtract UI and info panel
                    vp.canvasW -= INFO_W;
                    vp.canvasH -= UI_H;
                    vp.canvasY  = UI_H;
                }
                break;
            }
        }

        // ── Render ────────────────────────────────────────────────────────────

        // Background
        SDL_SetRenderDrawColor(renderer,
            COL_BG.r, COL_BG.g, COL_BG.b, 255);
        SDL_RenderClear(renderer);

        // Canvas clip rect
        SDL_Rect canvasRect{vp.canvasX, vp.canvasY, vp.canvasW, vp.canvasH};
        SDL_RenderSetClipRect(renderer, &canvasRect);

        // Draw triangulation
        drawTriangulation(renderer, state, vp);
        drawCorridor(renderer, state, vp);
        drawPoints(renderer, state, vp);
        drawPath(renderer, state, vp);

        // Start / goal markers
        if (state.startPt) {
            drawMarker(renderer, tc,
                       vp.toSX(state.startPt->x), vp.toSY(state.startPt->y),
                       COL_START, "S");
        }
        if (state.goalPt) {
            drawMarker(renderer, tc,
                       vp.toSX(state.goalPt->x), vp.toSY(state.goalPt->y),
                       COL_GOAL, "G");
        }

        SDL_RenderSetClipRect(renderer, nullptr);

        // ── UI panel ──────────────────────────────────────────────────────────
        fillRect(renderer, 0, 0, WIN_W, UI_H, COL_UI_BG);
        setColor(renderer, {170,180,200,255});
        SDL_RenderDrawLine(renderer, 0, UI_H-1, WIN_W, UI_H-1);

        int mx, my;
        SDL_GetMouseState(&mx, &my);

        filenameInput.draw(renderer, tc);
        btnLoad.draw(renderer, tc, btnLoad.contains(mx,my));
        btnDela.active = (state.algo == AlgoChoice::Delaunay);
        btnGrdy.active = (state.algo == AlgoChoice::Greedy);
        btnQgrd.active = (state.algo == AlgoChoice::Quasi);
        btnDela.draw(renderer, tc, btnDela.contains(mx,my));
        btnGrdy.draw(renderer, tc, btnGrdy.contains(mx,my));
        btnQgrd.draw(renderer, tc, btnQgrd.contains(mx,my));
        btnClear.draw(renderer, tc, btnClear.contains(mx,my));

        // Status message
        if (!state.statusMsg.empty()) {
            SDL_Color sc = state.statusIsError ? COL_ERROR : SDL_Color{30,100,30,255};
            tc.draw(state.statusMsg, 715, 18, sc);
        }

        // ── Info panel ────────────────────────────────────────────────────────
        drawInfoPanel(renderer, tc, state);

        // Canvas border
        setColor(renderer, {180,190,205,255});
        SDL_RenderDrawRect(renderer, &canvasRect);

        SDL_RenderPresent(renderer);
    }

    SDL_StopTextInput();
    if (font) TTF_CloseFont(font);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
